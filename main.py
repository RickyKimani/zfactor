from pathlib import Path
import json

from scripts import extractor
from scripts.parse import parse_char_prop, parse_lee_kesler_tables, parse_antoine_table


def main():
    pdf_source = "text.pdf"
    appendix_pdf = "data/appendix_BCD.pdf"
    b1_pdf = "data/b1_char_prop.pdf"
    b2_pdf = "data/b2_antoine.pdf"
    lee_kesler_pdf = "data/lee_kesler.pdf"
    json_out_lee_kesler = "data/lee_kesler.json"
    json_out_b1 = "data/b1_char_prop.json"
    json_out_b2 = "data/b2_antoine.json"

    extractor.extract_pages(
        input_pdf=pdf_source,
        start=682,
        end=711,
        out=appendix_pdf,
    )

    #properties
    extractor.extract_pages(
        input_pdf=appendix_pdf,
        start=1,
        end=3,
        out=b1_pdf
    )

    #antoine
    extractor.extract_pages(
        input_pdf=appendix_pdf,
        start=4,
        end=5,
        out=b2_pdf
    )

    #lee-kesler
    extractor.extract_pages(
        input_pdf=appendix_pdf,
        start=14,
        end=30,
        out=lee_kesler_pdf
    )

    data_properties = parse_char_prop(b1_pdf)
    Path(json_out_b1).parent.mkdir(parents=True, exist_ok=True)
    with open(json_out_b1, "w") as f:
        json.dump(data_properties, f, indent=2)
    print(f"Extracted {len(data_properties)} substances -> {json_out_b1}")
    data_antoine = parse_antoine_table(b2_pdf)
    with open(json_out_b2, "w") as f:
        json.dump(data_antoine, f, indent=2)
    print(f"Extracted {len(data_antoine)} antoine entries -> {json_out_b2}")
    data_lee = parse_lee_kesler_tables(lee_kesler_pdf)
    Path(json_out_lee_kesler).parent.mkdir(parents=True, exist_ok=True)
    with open(json_out_lee_kesler, "w") as f:
        json.dump(data_lee, f, indent=2)
    print(f"Extracted {len(data_lee)} tables -> {json_out_lee_kesler}")


if __name__ == "__main__":
    main()
