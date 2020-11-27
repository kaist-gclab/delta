# 용어 사전

## 에셋

- Asset: 에셋
- AssetTag: 에셋 태그
- AssetContent: 에셋 내용
  - Model
  - RenderedImage
  - Dictionary
- AssetType: 에셋 형식
  - text/plain
  - application/json
  - image/png
  - image/bmp
  - model/x.stl-ascii
  - model/x.stl-binary
  - model/obj
  - model/delta
- EncryptionKey: 암호화 키
  - Enabled: 활성 여부

## 뷰어

- Viewer: 뷰어

## 처리기

- ProcessorNode: 처리기 노드
  - node_a
  - node_b
  - ...
- ProcessorNodeStatus: 처리기 노드 상태
- ProcessorNodeCapability: 처리기 노드 능력
  - ConvertToDelta, model/x.stl-ascii
  - Render, model/delta
  - Explain, model/delta
  - Explain, image/png
  - GetSize, image/png
  - GetSize, image/bmp

## 작업

- Job: 작업
- JobExecution: 작업 실행
- JobExecutionStatus: 작업 실행 상태
- JobType: 작업 유형
  - Render
  - ConvertToDelta
  - Explain
  - GetSize
