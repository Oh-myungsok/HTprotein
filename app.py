import streamlit as st
from Bio import SeqIO
import pandas as pd
import os

# FASTA 파일 경로 (app.py와 같은 위치에 있어야 함)
FASTA_FILE = "d1_fasta_clean.fasta"

# 아미노산 분자량 대략값
aa_weights = {
    'A': 89.1, 'R': 174.2, 'N': 132.1, 'D': 133.1,
    'C': 121.2, 'E': 147.1, 'Q': 146.2, 'G': 75.1,
    'H': 155.2, 'I': 131.2, 'L': 131.2, 'K': 146.2,
    'M': 149.2, 'F': 165.2, 'P': 115.1, 'S': 105.1,
    'T': 119.1, 'W': 204.2, 'Y': 181.2, 'V': 117.1
}

def calc_mw(seq):
    return sum(aa_weights.get(aa, 0) for aa in seq)

# Streamlit UI
st.title("HotProtein Search App 🔬")
st.write("FASTA 파일에서 단백질을 검색합니다.")

# 사용자 입력
min_mw = st.number_input("최소 분자량", value=10000.0)
max_mw = st.number_input("최대 분자량", value=50000.0)
keyword = st.text_input("검색 키워드").lower()

if st.button("검색 실행"):
    if not os.path.exists(FASTA_FILE):
        st.error(f"FASTA 파일 {FASTA_FILE} 을 찾을 수 없습니다. 저장소 루트에 파일을 두세요.")
    else:
        results = []
        try:
            for record in SeqIO.parse(FASTA_FILE, "fasta"):
                mw = calc_mw(str(record.seq))
                if min_mw <= mw <= max_mw:
                    if keyword in record.description.lower():
                        results.append((record.id, mw, str(record.seq)[:50] + "..."))
        except Exception as e:
            st.error(f"FASTA 파일을 읽는 중 오류 발생: {e}")

        if results:
            st.subheader("검색 결과")
            st.write("총 결과 수:", len(results))
            df = pd.DataFrame(results, columns=["ID", "분자량", "서열(앞 50aa)"])
            st.dataframe(df)
        else:
            st.warning("조건에 맞는 단백질이 없습니다.")
