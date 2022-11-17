# 용어 사전

## 애셋

- Asset: 애셋
- AssetTag: 애셋 태그
- AssetType: 애셋 유형
  - Model
  - RenderedImage
  - Image
  - Dictionary
- MediaType: 미디어 형식
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
  - AssetContent
  - MediaType

## 처리기

- ProcessorNode: 처리기 노드
  - node_a
  - node_b
  - ...
- ProcessorNodeStatus: 처리기 노드 상태
- ProcessorNodeCapability: 처리기 노드 능력
  - ConvertToDelta, Model, model/x.stl-ascii
  - Render, Model, model/delta
  - Explain, Model, model/delta
  - GetSize, null, image/png
  - GetSize, null, image/bmp

## 작업

- Job: 작업
- JobExecution: 작업 실행
- JobExecutionStatus: 작업 실행 상태
- JobType: 작업 유형
  - Render
  - ConvertToDelta
  - Explain
  - GetSize
