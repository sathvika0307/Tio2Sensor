import markdown
from weasyprint import HTML, CSS

# Read the Markdown file
with open('final_report.md', 'r') as f:
    md_content = f.read()

# Convert Markdown to HTML
html_content = markdown.markdown(md_content, extensions=['tables'])

# Add some CSS for better formatting
css = CSS(string='''
body { font-family: Arial, sans-serif; margin: 40px; }
h1, h2, h3 { color: #333; }
table { border-collapse: collapse; width: 100%; }
th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
th { background-color: #f2f2f2; }
''')

# Convert HTML to PDF
HTML(string=html_content).write_pdf('final_report.pdf', stylesheets=[css])

print("Converted final_report.md to final_report.pdf")
