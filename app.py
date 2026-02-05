import streamlit as st
from Bio import SeqIO
import pandas as pd
import os
from io import BytesIO
import pickle
import glob

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

# 🔥 모든 classification.pkl 파일에서 train/test names & labels 읽어서 매핑 딕셔너리 생성
classification = {}
for file in glob.glob("d1_*_classification.pkl"):
    with open(file, "rb") as f:
        data = pickle.load(f)
        if "train_names" in data and "train_labels" in data:
            classification.update(dict(zip(data["train_names"], data["train_labels"])))
        if "test_names" in data and "test_labels" in data:
            classification.update(dict(zip(data["test_names"], data["test_labels"])))

# 🔎 FASTA ID 전처리 함수
def normalize_id(fasta_id):
    if "|" in fasta_id:
        parts = fasta_id.split("|")
        if len(parts) > 1:
            return parts[1]
    return fasta_id

# 🔎 클래스별 온도 범위 정의
temp_ranges = {
    0: "≤ 15°C",
    1: "≤ 20°C",
    2: "20–45°C",
    3: "45–60°C",
    4: "≥ 70°C",
    5: "≥ 80°C",
    6: "90–100°C 이상",
    7: "연구자 정의",
    8: "N/A",
    9: "N/A"
}

st.title("HotProtein Search App 🔬")
st.write("FASTA 파일에서 단백질을 검색합니다.")

# 🔎 검색 조건 입력 (기본값 변경)
min_mw = st.number_input("최소 분자량 (Da)", value=3000.0)   # 3KDa
max_mw = st.number_input("최대 분자량 (Da)", value=120000.0) # 120KDa
keyword = st.text_input("검색 키워드").lower()

# 🔎 열안정성 클래스 선택
selected_class = st.selectbox(
    "열안정성 클래스 선택 (8,9 제외)",
    options=[None,0,1,2,3,4,5,6,7],
    format_func=lambda x: "전체" if x is None else f"{x} 클래스"
)

# 📌 열안정성 클래스 설명 안내문 (표 형식)
st.markdown("""
**열안정성 클래스 설명 (0~9):**

| 클래스 | 구분 | 안정 온도 범위 | 설명 |
|--------|------|----------------|------|
| 0 | 극저온성 (Psychrophilic) | ≤ 15°C | 극저온 환경에서 성장 |
| 1 | 저온성 (Low-temperature) | ≤ 20°C | 저온 환경에서 안정 |
| 2 | 중온성 (Mesophilic) | 20–45°C | 대부분 생물 성장 온도 |
| 3 | 약간 고온성 (Moderately thermophilic) | 45–60°C | 온천 등 고온 환경 |
| 4 | 고온성 (Thermophilic) | ≥ 70°C | 고온 환경에서 안정 |
| 5 | 초고온성 (Hyperthermophilic) | ≥ 80°C | 심해 열수구 등 극한 고온 |
| 6 | 극한 환경성 (Extreme environment) | 90–100°C 이상 | 극한 조건에서 안정 |
| 7 | 변형된 안정성 (Engineered stability) | 연구자 정의 | 인위적 변이로 안정성 강화 |
| 8 | 불안정성 (Unstable) | N/A | 안정성 낮음 |
| 9 | 기타/분류 불명 (Miscellaneous) | N/A | 데이터 부족 |
""")

if st.button("검색 실행"):
    if not os.path.exists(FASTA_FILE):
        st.error(f"FASTA 파일 {FASTA_FILE} 을 찾을 수 없습니다.")
    else:
        results_display = []
        results_save = []

        try:
            for record in SeqIO.parse(FASTA_FILE, "fasta"):
                mw = calc_mw(str(record.seq))
                if min_mw <= mw <= max_mw:
                    if keyword in record.description.lower():
                        accession = normalize_id(record.id)
                        thermo_class = classification.get(accession, "N/A")
                        temp_range = temp_ranges.get(thermo_class, "N/A")

                        if selected_class is None or thermo_class == selected_class:
                            results_display.append((record.id, round(mw, 2), str(record.seq)[:100] + "...", thermo_class, temp_range))
                            results_save.append((record.id, round(mw, 2), str(record.seq), thermo_class, temp_range))
        except Exception as e:
            st.error(f"FASTA 파일을 읽는 중 오류 발생: {e}")

        if results_display:
            st.subheader("검색 결과 (앞 100aa 표시)")
            st.write("총 결과 수:", len(results_display))

            df_display = pd.DataFrame(results_display, columns=["ID", "분자량(Da)", "서열(앞 100aa)", "열안정성 클래스", "온도 범위"])
            df_save = pd.DataFrame(results_save, columns=["ID", "분자량(Da)", "서열 전체", "열안정성 클래스", "온도 범위"])

            def highlight_class(val):
                if val == 0:
                    return 'color: blue; font-weight: bold;'
                elif val in [4,5]:
                    return 'color: red; font-weight: bold;'
                elif val == 2:
                    return 'color: green; font-weight: bold;'
                else:
                    return ''
            
            styled_df = df_display.style.applymap(highlight_class, subset=["열안정성 클래스"])
            st.dataframe(styled_df, height=400, use_container_width=True)

            csv_data = df_save.to_csv(index=False).encode("utf-8")
            st.download_button("📥 CSV 파일로 저장", csv_data, "search_results.csv", "text/csv")

            buffer = BytesIO()
            with pd.ExcelWriter(buffer, engine="xlsxwriter") as writer:
                df_save.to_excel(writer, index=False, sheet_name="Results")
            st.download_button("📥 엑셀 파일로 저장", buffer.getvalue(),
                               "search_results.xlsx",
                               "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
        else:
            st.warning("조건에 맞는 단백질이 없습니다.")
