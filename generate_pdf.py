from fpdf import FPDF
import markdown

# Read the Markdown file
with open('final_report.md', 'r') as f:
    md_content = f.read()

# Convert Markdown to HTML
html_content = markdown.markdown(md_content)

# Create PDF
pdf = FPDF()
pdf.add_page()
pdf.set_font("Arial", size=12)

# Simple text addition (since HTML parsing is complex, just add raw text)
lines = md_content.split('\n')
for line in lines:
    if line.startswith('# '):
        pdf.set_font("Arial", 'B', 16)
        pdf.cell(200, 10, txt=line[2:], ln=True)
        pdf.set_font("Arial", size=12)
    elif line.startswith('## '):
        pdf.set_font("Arial", 'B', 14)
        pdf.cell(200, 10, txt=line[3:], ln=True)
        pdf.set_font("Arial", size=12)
    elif line.startswith('|'):
        pdf.cell(200, 10, txt=line, ln=True)
    elif line.strip():
        pdf.cell(200, 10, txt=line, ln=True)
    else:
        pdf.ln(5)

pdf.output("final_report.pdf")
print("PDF generated as final_report.pdf")
