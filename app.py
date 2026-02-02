import streamlit as st
from Bio import SeqIO
from io import StringIO
import plotly.express as px
import pandas as pd
import re
from fpdf import FPDF  # برای PDF

st.set_page_config(page_title="FASTA Analyzer Pro", layout="wide")
st.title("🚀 FASTA Analyzer Pro – بیوتک واقعی شروع شد")
st.markdown("آپلود چند فایل • جستجوی motif • دانلود گزارش CSV/PDF")

uploaded_files = st.file_uploader(
    "فایل(های) FASTA رو آپلود کن",
    type=["fasta", "fa", "fas", "fna", "fasta.gz", "txt"],
    accept_multiple_files=True
)

motif_pattern = st.text_input(
    "جستجوی Motif (مثل ATG یا regex ساده مثل A{3,} یا [ATGC]{5}):",
    placeholder="ATG"
)

if uploaded_files:
    all_records = []
    for uploaded in uploaded_files:
        try:
            string_data = StringIO(uploaded.getvalue().decode("utf-8"))
            records = list(SeqIO.parse(string_data, "fasta"))
            all_records.extend(records)
            st.info(f"✓ {uploaded.name} → {len(records)} سکانس")
        except Exception as e:
            st.error(f"خطا در خواندن {uploaded.name}: {e}")

    if all_records:
        st.success(f"کلاً {len(all_records)} سکانس آماده! ")

        # پیام وضعیت جستجوی motif
        if motif_pattern:
            st.info(f"جستجو برای motif: **{motif_pattern.upper()}** (case-insensitive)")
        else:
            st.info("هیچ motifi وارد نشده — فقط آمار پایه نمایش داده می‌شود")

        # آماده کردن داده برای جدول و گزارش
        data = []
        for i, rec in enumerate(all_records):
            seq = str(rec.seq).upper()
            length = len(seq)
            gc = (seq.count('G') + seq.count('C')) / length * 100 if length > 0 else 0
            counts = {
                'A': seq.count('A'),
                'T': seq.count('T'),
                'G': seq.count('G'),
                'C': seq.count('C'),
                'N': seq.count('N')
            }

            # جستجوی motif — case-insensitive و 1-based
            motif_info = "پیدا نشد"
            positions = []
            if motif_pattern:
                try:
                    # re.compile برای regex کامل و case-insensitive
                    pattern_compiled = re.compile(motif_pattern, re.IGNORECASE)
                    matches = list(pattern_compiled.finditer(seq))
                    positions = [m.start() + 1 for m in matches]  # 1-based
                    if positions:
                        matched_seqs = [m.group() for m in matches]
                        motif_info = ", ".join(f"{pos} ({matched_seq})" for pos, matched_seq in zip(positions, matched_seqs))
                except re.error:
                    motif_info = "خطا در regex — pattern نامعتبر"

            data.append({
                'ID': rec.id,
                'Length': length,
                'GC%': round(gc, 2),
                'A': counts['A'],
                'T': counts['T'],
                'G': counts['G'],
                'C': counts['C'],
                'N': counts['N'],
                'Motif Matches': motif_info
            })

            # نمایش جزئیات هر سکانس در expander
            with st.expander(f"سکانس {i+1}: {rec.id} – طول: {length:,} nt", expanded=(i < 3)):
                if length > 500:
                    st.code(str(rec.seq)[:500] + "...")
                else:
                    st.code(str(rec.seq))

                col1, col2 = st.columns(2)
                col1.metric("طول", f"{length:,} nt")
                col2.metric("GC%", f"{gc:.2f}%")

                df_counts = pd.DataFrame.from_dict(counts, orient='index', columns=['تعداد'])
                fig = px.bar(
                    df_counts,
                    text_auto=True,
                    color=df_counts.index,
                    color_discrete_sequence=px.colors.qualitative.Bold
                )
                fig.update_layout(title="توزیع نوکلئوتید", showlegend=False, height=400)
                st.plotly_chart(fig, use_container_width=True)

                if motif_pattern:
                    st.write(f"موتيف '{motif_pattern.upper()}': {motif_info}")

        # جدول خلاصه همه سکانس‌ها
        df = pd.DataFrame(data)
        st.subheader("📊 جدول خلاصه همه سکانس‌ها")
        st.dataframe(
            df.style.format({
                "Length": "{:,}",
                "GC%": "{:.2f}"
            }),
            use_container_width=True
        )

        # دانلود CSV
        csv = df.to_csv(index=False).encode('utf-8')
        st.download_button(
            label="📥 دانلود گزارش کامل به صورت CSV",
            data=csv,
            file_name="FASTA_analysis_report.csv",
            mime="text/csv"
        )

        # دانلود PDF
        def generate_pdf(report_df):
            pdf = FPDF()
            pdf.add_page()
            pdf.set_font("Helvetica", size=16)
            pdf.cell(200, 10, txt="FASTA Analyzer Pro Report", ln=1, align='C')
            pdf.ln(10)
            pdf.set_font("Helvetica", size=10)

            # هدر جدول
            col_widths = [40, 20, 20, 15, 15, 15, 15, 15, 50]
            headers = list(report_df.columns)
            for j, header in enumerate(headers):
                pdf.cell(col_widths[j], 10, header, 1)
            pdf.ln()

            # داده‌ها
            for _, row in report_df.iterrows():
                for j, val in enumerate(row):
                    text = str(val)
                    if len(text) > 20:
                        text = text[:17] + "..."
                    pdf.cell(col_widths[j], 10, text, 1)
                pdf.ln()

            return pdf.output(dest="S").encode("latin-1")

        if st.button("📄 تولید و دانلود گزارش PDF"):
            with st.spinner("در حال ساخت PDF..."):
                pdf_data = generate_pdf(df)
            st.download_button(
                label="📥 دانلود PDF حالا",
                data=pdf_data,
                file_name="FASTA_analysis_report.pdf",
                mime="application/pdf"
            )

else:
    st.info("هنوز فایلی آپلود نشده. منتظر فایل‌های FASTA شما هستیم! ")
