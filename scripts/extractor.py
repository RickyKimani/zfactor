from pathlib import Path
from PyPDF2 import PdfReader, PdfWriter


def extract_pages(
    input_pdf: str,
    start: int,
    end: int,
    out: str = "data/extracted.pdf",
) -> None:
    """
    Extract a range of pages from a PDF file and save them to a new PDF.

    This function takes an input PDF, extracts pages from 'start' to 'end' (inclusive),
    and writes them to a new PDF file. Page numbers are 1-based (e.g., start=1 for the first page).

    Parameters:
    - input_pdf (str): Path to the input PDF file.
    - start (int): Starting page number (1-based).
    - end (int): Ending page number (1-based).
    - out (str): Path to the output PDF file. Defaults to "data/extracted.pdf".

    Raises:
    - FileNotFoundError: If the input PDF does not exist.
    - ValueError: If start or end are out of range for the PDF.

    Example:
    extract_pages("book.pdf", 696, 711, "data/lee_kesler_tables.pdf")
    """
    input_path = Path(input_pdf)
    output_path = Path(out)

    if not input_path.exists():
        raise FileNotFoundError(f"Input PDF not found: {input_path}")

    start_idx = start - 1
    end_idx = end - 1

    print(f"Extracting pages {start} - {end} from {input_pdf}")

    reader = PdfReader(input_path)
    writer = PdfWriter()

    num_pages = len(reader.pages)

    if start_idx >= num_pages or end_idx >= num_pages:
        raise ValueError(f"PDF has only {num_pages} pages!")

    for i in range(start_idx, end_idx + 1):
        writer.add_page(reader.pages[i])

    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "wb") as f:
        writer.write(f)

    print(f"Saved extracted pages to {output_path.resolve()}")

if __name__ == "__main__":
    extract_pages("./text.pdf", 682, 711)