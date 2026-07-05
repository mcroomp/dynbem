import fitz  # PyMuPDF
import sys

pdf_path = r"C:\Users\krist\Downloads\ADA099532.pdf"
output_path = r"C:\Users\krist\Downloads\ADA099532.txt"

try:
    doc = fitz.open(pdf_path)
    text_output = []
    
    for page_num in range(len(doc)):
        page = doc[page_num]
        text = page.get_text()
        text_output.append(f"--- PAGE {page_num + 1} ---\n{text}\n")
    
    full_text = "\n".join(text_output)
    
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(full_text)
    
    print(f"✓ Extracted {len(doc)} pages")
    print(f"✓ Output: {output_path}")
    print(f"✓ Size: {len(full_text)} characters")
    
except Exception as e:
    print(f"Error: {e}", file=sys.stderr)
    sys.exit(1)
