import streamlit as st
from Bio import SeqIO
import pandas as pd
import os
from io import BytesIO

FASTA_FILE = "d1_fasta_clean.fasta"

aa_weights = {
    'A': 89.1, 'R': 174.2, 'N': 132.1, 'D': 133.1,
    'C': 121.2, 'E': 147.1, 'Q': 146.2, 'G': 75.1,
    'H': 155.2, 'I': 131.2, 'L': 131.2, 'K': 146.2,
    'M': 149.2, 'F': 165.2, 'P': 115.1, 'S': 105.1,
    'T': 119.1, 'W': 204.2, 'Y': 181.2, 'V': 117.1
}

def calc_mw(seq):
    total = sum(aa_weights.get(aa, 0) for aa in seq)
    water_loss = (len(seq) - 1) * 18.0 if len(seq) > 1 else 0
    return total - water_loss

st.title("HotProtein Search App 🔬")
st.write("FASTA 파일에서 단백질을 검색합니다.")

min_mw = st.number_input("최소 분자량 (Da)", value=10000.0)
max_mw = st.number_input("최대 분자량 (Da)", value=50000.0)
keyword = st.text_input("검색 키워드").lower()

if st.button("검색 실행"):
    if not os.path.exists(FASTA_FILE):
        st.error(f"FASTA 파일 {FASTA_FILE} 을 찾을 수 없습니다.")
    else:
        results_display = []  # 화면 표시용 (앞 100aa)
        results_save = []     # 저장용 (전체 서열)

        try:
            for record in SeqIO.parse(FASTA_FILE, "fasta"):
                mw = calc_mw(str(record.seq))
                if min_mw <= mw <= max_mw:
                    if keyword in record.description.lower():
                        # 화면 표시용: 앞 100개만
                        results_display.append((record.id, round(mw, 2), str(record.seq)[:100] + "..."))
                        # 저장용: 전체 서열
                        results_save.append((record.id, round(mw, 2), str(record.seq)))
        except Exception as e:
            st.error(f"FASTA 파일을 읽는 중 오류 발생: {e}")

        if results_display:
            st.subheader("검색 결과 (앞 100aa 표시)")
            st.write("총 결과 수:", len(results_display))

            df_display = pd.DataFrame(results_display, columns=["ID", "분자량(Da)", "서열(앞 100aa)"])
            df_save = pd.DataFrame(results_save, columns=["ID", "분자량(Da)", "서열 전체"])

            # 📌 세로 스크롤만 가능하게, 가로 스크롤 제거
            st.dataframe(df_display, height=400, use_container_width=True)

            # 📥 CSV 다운로드 버튼 (전체 서열)
            csv_data = df_save.to_csv(index=False).encode("utf-8")
            st.download_button(
                label="📥 CSV 파일로 저장 (전체 서열)",
                data=csv_data,
                file_name="search_results.csv",
                mime="text/csv"
            )

            # 📥 엑셀 다운로드 버튼 (전체 서열)
            buffer = BytesIO()
            with pd.ExcelWriter(buffer, engine="xlsxwriter") as writer:
                df_save.to_excel(writer, index=False, sheet_name="Results")
            st.download_button(
                label="📥 엑셀 파일로 저장 (전체 서열)",
                data=buffer.getvalue(),
                file_name="search_results.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            )

        else:
            st.warning("조건에 맞는 단백질이 없습니다.")
