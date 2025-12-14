import streamlit as st
from Bio import SeqIO
import plotly.express as px
import pandas as pd

st.set_page_config(page_title="FASTA Analyzer", layout="wide")
st.title("FASTA Analyzer – بیوتک واقعی شروع شد")

uploaded = st.file_uploader("فایل FASTA خودت رو اینجا آپلود کن", type=["fasta", "fa", "fas"])

if uploaded:
    records = list(SeqIO.parse(uploaded, "fasta"))
    st.success(f"{len(records)} سکانس پیدا شد!")
    
    for i, rec in enumerate(records):
        st.subheader(f"سکانس {i+1}: {rec.id}")
        st.code(str(rec.seq)[:500] + "..." if len(rec.seq) > 500 else rec.seq)
        
        seq = str(rec.seq)
        length = len(seq)
        gc = (seq.count('G') + seq.count('C')) / length * 100
        
        st.metric("طول سکانس", length)
        st.metric("GC Content", f"{gc:.2f}%")
        
        counts = {'A': seq.count('A'), 'T': seq.count('T'), 
                  'G': seq.count('G'), 'C': seq.count('C'), 'N': seq.count('N')}
        df = pd.DataFrame.from_dict(counts, orient='index', columns=['تعداد'])
        fig = px.bar(df, text_auto=True, title="توزیع نوکلئوتید")
        st.plotly_chart(fig, use_container_width=True)

st.write("نسخه ۱ — به زودی جستجوی motif و گزارش PDF هم میاد")/ تک به تک بگو اینو چه جوری انجامش بدم هیچی نمیفهمم ازش
