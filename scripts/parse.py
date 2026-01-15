from pathlib import Path
from PyPDF2 import PdfReader
import re

FLOAT = r"[-+]?(?:\d*\.\d+|\d+)"


def is_data_row(line: str) -> bool:
    """
    A valid data row:
    - starts with a letter or digit (substance name)
    - contains at least 6 numeric fields
    """
    nums = re.findall(FLOAT, line)
    return len(nums) >= 6 and (line[0].isalpha() or line[0].isdigit())


def parse_char_prop(pdf_path: str) -> list:
    """
    Parse Table B.1 (Characteristics of Pure Substances).

    Assumes:
    - PDF contains ONLY the table pages
    - No appendix boundary detection needed
    """
    reader = PdfReader(pdf_path)
    substances = []

    for page in reader.pages:
        text = page.extract_text()
        if not text:
            continue

        lines = [l.strip() for l in text.splitlines() if l.strip()]

        for line in lines:
            # Clean up special characters
            line = line.replace("†", "").replace("*", "")
            # Normalize unicode minus signs to ASCII hyphen
            line = line.replace("\u2212", "-").replace("\u2013", "-")

            # Skip known headers / footers
            if (
                "Table B.1" in line
                or "Molar mass" in line
                or "Critical properties" in line
                or "Substance" in line
            ):
                continue

            if not is_data_row(line):
                continue

            matches = list(re.finditer(FLOAT, line))
            if len(matches) < 6:
                continue

            # --- Special case: 6 columns (Missing Omega OR Missing Tn) ---
            if len(matches) == 6:
                name_end = matches[0].start()
                name = line[:name_end].strip()

                nums = [m.group() for m in matches]

                try:
                    vals = [float(n) for n in nums]

                    # Heuristic: Check 2nd number.
                    # If > 10, it is likely Tc (Missing Omega)
                    # If <= 10, it is likely Omega (Missing Tn)
                    if vals[1] > 10.0:
                        mw = vals[0]
                        omega = 0.0
                        tc = vals[1]
                        pc = vals[2]
                        zc = vals[3]
                        vc = vals[4]
                        tn = vals[5]
                    else:
                        mw = vals[0]
                        omega = vals[1]
                        tc = vals[2]
                        pc = vals[3]
                        zc = vals[4]
                        vc = vals[5]
                        tn = 0.0
                except ValueError:
                    continue

            # --- Normal case: full 7-column row ---
            else:
                data_matches = matches[-7:]
                name_end = data_matches[0].start()
                name = line[:name_end].strip()

                nums = [m.group() for m in data_matches]

                try:
                    mw = float(nums[0])
                    omega = float(nums[1])
                    tc = float(nums[2])
                    pc = float(nums[3])
                    zc = float(nums[4])
                    vc = float(nums[5])
                    tn = float(nums[6])
                except ValueError:
                    continue

            substances.append({
                "name": name,
                "mw": mw,
                "acentric": omega,
                "tn": tn,
                "critical": {
                    "tc": tc,
                    "pc": pc,
                    "vc": vc,
                    "zc": zc,
                },
            })

    return substances

def parse_lee_kesler_tables(pdf_path: str) -> dict:
    path = Path(pdf_path)
    if not path.exists():
        raise FileNotFoundError(f"{pdf_path} not found")

    reader = PdfReader(path)
    all_tables = {}

    table_keys = {
        "Z0": "z0",
        "Z 1": "z1",
        "Z1": "z1",
        "( HR)0": "h0",
        "( HR)1": "h1",
        "( SR)0": "s0",
        "( SR)1": "s1",
        "ϕ0": "phi0",
        "ϕ1": "phi1",
    }

    for page_num, page in enumerate(reader.pages):
        text = page.extract_text()
        lines = [l.strip() for l in text.splitlines() if l.strip()]

        table_title = None
        for line in lines:
            if "Table D." in line and ":" in line:
                table_title = line
                break

        if not table_title:
            continue

        key = None
        for k, v in table_keys.items():
            if k in table_title:
                key = v
                break

        if not key:
            continue

        pressures = []
        data_rows = []

        for line in lines:
            # Normalize unicode minus signs to ASCII hyphen
            line = line.replace("\u2212", "-").replace("\u2013", "-")

            if "Pr =" in line:
                pressures = [float(x) for x in re.findall(FLOAT, line)]
                continue

            nums = re.findall(FLOAT, line)
            if not nums or len(nums) < 2:
                continue

            try:
                tr = float(nums[0])
            except ValueError:
                continue

            z_vals = []
            for x in nums[1 : len(pressures) + 1]:
                try:
                    z_vals.append(float(x))
                except ValueError:
                    z_vals.append(None)

            if z_vals and len(z_vals) == len(pressures):
                data_rows.append((tr, z_vals))

        if pressures and data_rows:
            reduced_temperatures = [r[0] for r in data_rows]
            matrix = [r[1] for r in data_rows]

            table_data = {
                "reduced_pressure": pressures,
                "reduced_temperature": reduced_temperatures,
                "values": matrix,
            }

            if key not in all_tables:
                all_tables[key] = []
            all_tables[key].append(table_data)

            print(f"Parsed {table_title}: {len(data_rows)} rows")

    # Combine low and high Pr tables
    for key in all_tables:
        if len(all_tables[key]) == 2:
            t1, t2 = all_tables[key]
            pr1 = t1["reduced_pressure"]
            pr2 = t2["reduced_pressure"]
            combined_pr = pr1 + pr2[1:]
            combined_values = []
            for r in range(len(t1["values"])):
                row = t1["values"][r] + t2["values"][r][1:]
                combined_values.append(row)
            combined_table = {
                "reduced_pressure": combined_pr,
                "reduced_temperature": t1["reduced_temperature"],
                "values": combined_values,
            }
            all_tables[key] = [combined_table]

    return all_tables

def parse_antoine_table(pdf_path: str) -> list:
    reader = PdfReader(pdf_path)
    substances = []

    for page in reader.pages:
        text = page.extract_text()
        if not text:
            continue

        lines = [l.strip() for l in text.splitlines() if l.strip()]

        for line in lines:
            # Skip headers
            if "Table B.2" in line or "Constants" in line or "Formula" in line:
                continue

            # Clean up special characters
            line = line.replace("†", "").replace("*", "")
            # Normalize unicode minus signs and separators
            line = line.replace("\u2212", "-").replace("\u2013", "-").replace("\u2014", " ")

            # Find numbers
            # We expect 7 numbers: A, B, C, tmin, tmax, H, tn
            nums = re.findall(FLOAT, line)
            if len(nums) < 7:
                continue

            # The last 7 matches are likely our data
            data_matches = list(re.finditer(FLOAT, line))[-7:]

            # Name and Formula are before the first number (A)
            first_num_idx = data_matches[0].start()
            prefix = line[:first_num_idx].strip()

            # Split prefix into Name and Formula
            # Formula is usually the last part of the prefix
            parts = prefix.split()
            if len(parts) < 2:
                continue

            # Fix for detached formula subscripts (e.g. "CCl 4" -> "CCl4")
            if len(parts) > 2 and parts[-1].isdigit():
                parts[-2] = parts[-2] + parts[-1]
                parts.pop()

            formula = parts[-1]
            name = " ".join(parts[:-1])

            vals = [float(m.group()) for m in data_matches]

            substances.append({
                "name": name,
                "formula": formula,
                "a": vals[0],
                "b": vals[1],
                "c": vals[2],
                "t_min": vals[3],
                "t_max": vals[4],
                "h": vals[5],
                "tn": vals[6]
            })

    return substances