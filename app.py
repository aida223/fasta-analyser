import streamlit as st
from Bio import SeqIO
from io import StringIO  # این خط جدید اضافه شده
import plotly.express as px
import pandas as pd

st.set_page_config(page_title="FASTA Analyzer", layout="wide")
st.title("🚀 FASTA Analyzer – بیوتک واقعی شروع شد")
st.markdown("چند فایل FASTA همزمان آپلود کن • تا ۲۰۰ مگابایت • کاملاً فارسی")

uploaded_files = st.file_uploader(
    "فایل(های) FASTA رو اینجا بکش یا انتخاب کن (چندتا هم می‌تونی همزمان)",
    type=["fasta", "fa", "fas", "fna", "fasta.gz", "txt"],
    accept_multiple_files=True
)

if uploaded_files:
    all_records = []
    for uploaded in uploaded_files:
        try:
            # اینجا تبدیل به text می‌کنیم (مهم‌ترین تغییر!)
            string_data = StringIO(uploaded.getvalue().decode("utf-8"))
            records = list(SeqIO.parse(string_data, "fasta"))
            all_records.extend(records)
            st.info(f"✓ {uploaded.name} → {len(records)} سکانس خونده شد")
        except Exception as e:
            st.error(f"خطا در خوندن {uploaded.name}: {e}")

    if all_records:
        st.success(f"کلاً {len(all_records)} سکانس آماده تحلیل شد! 🎉")

        for i, rec in enumerate(all_records):
            with st.expander(f"سکانس {i+1}: {rec.id} – طول: {len(rec.seq):,} nt", expanded=(i < 3)):
                if len(rec.seq) > 500:
                    st.code(str(rec.seq)[:500] + "... (ادامه سکانس مخفی شد)")
                else:
                    st.code(str(rec.seq))

                seq = str(rec.seq).upper()
                length = len(seq)
                gc = (seq.count('G') + seq.count('C')) / length * 100 if length > 0 else 0

                col1, col2 = st.columns(2)
                col1.metric("طول سکانس", f"{length:,} nt")
                col2.metric("GC Content", f"{gc:.2f}%")

                counts = {'A': seq.count('A'), 'T': seq.count('T'), 
                          'G': seq.count('G'), 'C': seq.count('C'), 'N': seq.count('N')}
                df = pd.DataFrame.from_dict(counts, orient='index', columns=['تعداد'])
                fig = px.bar(df, text_auto=True, color=counts.keys(), 
                             color_discrete_sequence=px.colors.qualitative.Bold)
                fig.update_layout(title="توزیع نوکلئوتید", showlegend=False)
                st.plotly_chart(fig, use_container_width=True)


