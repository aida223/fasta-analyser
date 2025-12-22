import streamlit as st
from Bio import SeqIO
from io import StringIO
import plotly.express as px
import pandas as pd
import re  # برای جستجوی motif
from fpdf2 import FPDF  # برای PDF (باید به requirements.txt اضافه کنی)

st.set_page_config(page_title="FASTA Analyzer Pro", layout="wide")
st.title("🚀 FASTA Analyzer Pro – بیوتک واقعی شروع شد")
st.markdown("آپلود چند فایل • جستجوی motif • دانلود گزارش CSV/PDF")

uploaded_files = st.file_uploader(
    "فایل(های) FASTA رو آپلود کن",
    type=["fasta", "fa", "fas", "fna", "fasta.gz", "txt"],
    accept_multiple_files=True
)

motif_pattern = st.text_input("جستجوی Motif (مثل ATG یا regex ساده):", placeholder="ATG")

if uploaded_files:
    all_records = []
    for uploaded in uploaded_files:
        try:
            string_data = StringIO(uploaded.getvalue().decode("utf-8"))
            records = list(SeqIO.parse(string_data, "fasta"))
            all_records.extend(records)
            st.info(f"✓ {uploaded.name} → {len(records)} سکانس")
        except Exception as e:
            st.error(f"خطا: {e}")

    if all_records:
        st.success(f"کلاً {len(all_records)} سکانس آماده! 🎉")

        # آماده کردن داده برای CSV
        data = []
        for i, rec in enumerate(all_records):
            seq = str(rec.seq).upper()
            length = len(seq)
            gc = (seq.count('G') + seq.count('C')) / length * 100 if length > 0 else 0
            counts = {'A': seq.count('A'), 'T': seq.count('T'), 'G': seq.count('G'), 'C': seq.count('C'), 'N': seq.count('N')}
            
            # جستجوی motif
            positions = [m.start() for m in re.finditer(motif_pattern, seq)] if motif_pattern else []
            motif_info = f"موقعیت‌ها: {positions}" if positions else "پیدا نشد"

            data.append({
                'ID': rec.id, 'Length': length, 'GC%': gc, 
                'A': counts['A'], 'T': counts['T'], 'G': counts['G'], 'C': counts['C'], 'N': counts['N'],
                'Motif Positions': motif_info
            })

            with st.expander(f"سکانس {i+1}: {rec.id} – طول: {length:,} nt", expanded=(i < 3)):
                st.code(str(rec.seq)[:500] + "..." if length > 500 else str(rec.seq))
                col1, col2 = st.columns(2)
                col1.metric("طول", f"{length:,} nt")
                col2.metric("GC%", f"{gc:.2f}%")
                df_counts = pd.DataFrame.from_dict(counts, orient='index', columns=['تعداد'])
                fig = px.bar(df_counts, text_auto=True, color=counts.keys(), color_discrete_sequence=px.colors.qualitative.Bold)
                fig.update_layout(title="توزیع نوکلئوتید", showlegend=False)
                st.plotly_chart(fig, use_container_width=True)
                
                if motif_pattern:
                    st.write(f"موتيف '{motif_pattern}': {motif_info}")

        # دانلود CSV
        df = pd.DataFrame(data)
        csv = df.to_csv(index=False).encode('utf-8')
        st.download_button("دانلود گزارش CSV", csv, "fasta_report.csv", "text/csv")

        # دانلود PDF
        def create_pdf():
            pdf = FPDF()
            pdf.add_page()
            pdf.set_font("Arial", size=12)
            pdf.cell(200, 10, txt="FASTA Analyzer Report", ln=1, align='C')
            for row in data:
                pdf.cell(200, 10, txt=f"ID: {row['ID']} | Length: {row['Length']} | GC: {row['GC%']:.2f}%", ln=1)
                pdf.cell(200, 10, txt=f"Motif: {row['Motif Positions']}", ln=1)
            return pdf.output(dest="S").encode('latin1')
        
        pdf_bytes = create_pdf()
        st.download_button("دانلود گزارش PDF", pdf_bytes, "fasta_report.pdf", "application/pdf")



