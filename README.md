# Jupyter Notebook Collection

이 디렉터리는 여러 연구/분석용 Jupyter 노트북을 포함하고 있습니다. GitHub에 정리하기 위해 다음과 같은 기본 파일을 추가했습니다.

- `.gitignore`: `.ipynb_checkpoints/`와 임시 파일을 제외합니다.
- `README.md`: 현재 디렉터리의 노트북 목록과 구성 정보를 제공합니다.

## 폴더 구성

### Annotation ML
- `GC.ipynb`
- `GC mac final.ipynb`
- `GC ML mac.ipynb`
- `HumanMacrophage_select.ipynb`
- `Human_ML-test.ipynb`
- `Human_ML-test-Copy1.ipynb`
- `Machinelearning.ipynb`
- `Untitled.ipynb`

### CRC - 1
- `bulk spot.ipynb`
- `Cancergenome.ipynb`
- `CRC singlecell.ipynb`
- `DEC_DATAMERGE.ipynb`
- `deconvolution.ipynb`
- `Public.ipynb`
- `Untitled-Copy1.ipynb`
- `Untitled.ipynb`

### CRC - 2
- `miRNA.ipynb`
- `Myeloid Proportion.ipynb`
- `Myeloid subtyping.ipynb`
- `Untitled.ipynb`

### CRC(Spatial)
- `VisiumHD_GSM8594567_P1CRC.ipynb`
- `VisiumHD_GSM8594568_P2CRC.ipynb`
- `VisiumHD_GSM8594569_P5CRC.ipynb`
- `Visium_merge.ipynb`
- `data/Single_cell.ipynb`

### GC
- `GC datamerge.ipynb`
- `pyscnic.ipynb`

### Oval - Human
- `Prework.ipynb`
- `scRNA (subset).ipynb`
- `CCRCC/ccRCC.ipynb`
- `CCRCC/ccRCC (ST).ipynb`
- `CCRCC/ccRCC Tumor spot.ipynb`

### Oval - Mouse
- `Cell oracle.ipynb`
- `Pseudotime Robustness.ipynb`
- `RNA-seq/class.ipynb`
- `Single_cell.ipynb`
- `Data calling/Datacalling.ipynb`
- `Data calling/Single_cell.ipynb`
- `spatial/Deconvolution.ipynb`
- `spatial/Tgfb_pathway.ipynb`
- `spatial/Visium(Old&Young).ipynb`
- `spatial/Visium_total_expression.ipynb`
- `Uterus_Epithelial.ipynb`
- `Pseudotime perterbation test/Pseudotime_simulation.ipynb`

## GitHub에 올리는 방법

1. 이 디렉터리에서 `git init` 실행
2. `git add .`로 추적할 파일 추가
3. `git commit -m "Initial notebook organization"`
4. GitHub에 새 리포지터리를 만든 뒤 원격 주소를 연결하고 `git push` 실행

```bash
git init
git add .
git commit -m "Initial notebook organization"
git remote add origin <your-github-repo-url>
git branch -M main
git push -u origin main
```

## 다음 단계

- 원하시면 이 폴더를 Git 저장소로 초기화해드릴 수 있습니다.
- 노트북을 Python 스크립트(`.py`)로 변환하거나, 주요 분석별로 `README`를 더 상세히 만드는 작업도 가능합니다.
