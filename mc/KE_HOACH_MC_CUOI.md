# Kế hoạch MC cuối cho bản nộp (bản nháp để duyệt, 28/09/2026)

Nguyên tắc: **một phiên bản lệnh (đóng băng), một lần chạy, một bộ log.** Mọi
bảng trong bài lấy từ lần chạy này. Các lần chạy trước (Study A 0.9.35, Study B
0.9.36, supplement kink, Supplement 2/2b) là giai đoạn kiểm tra lệnh; chúng
không đưa vào gói nộp.

## 0. Việc phải xong trước khi đóng băng (giai đoạn kiểm tra)

| # | Việc | Ai | Trạng thái |
|---|---|---|---|
| 0.1 | Kết quả 2b: power ở cấu hình Gong–Seo có gần Bảng 2 của họ không? Nếu không, cân nhắc sửa bootstrap trước khi đóng băng | VPS + tôi | đang chạy |
| 0.2 | `td`: giữ trong help hay ẩn (đề xuất: giữ) | bạn | chờ |
| 0.3 | Sửa thông báo còn nhắc `static`/`conttest`; `e(continuity_test)` | tôi | chờ 0.2 |
| 0.4 | Mặc định của lệnh (`grid(100) gridci(100) boot(299) trim(.10) refine(0)`): giữ hay đổi | bạn + tôi | chờ |
| 0.5 | Rà code lần cuối + chạy lại test Stata trên bản cuối | tôi + VPS | chờ |
| 0.6 | Số phiên bản nộp: 0.9.37 hay 1.0.0; gắn tag git | bạn | chờ |
| 0.7 | Nâng harness: một hợp đồng phiên bản (bản đóng băng), bỏ cột continuity, thêm các khối mới, smoke toàn bộ | tôi + VPS | sau 0.1–0.6 |

## 1. Thiết lập chung

- DGP: benchmark Gong–Seo (2026) như bản thảo Mục 4.1; T = 6 trừ khi ghi khác.
- Hai cấu hình công cụ, chạy song song:
  - **L3** = `maxlag(1 3)` (giới hạn độ trễ; cấu hình của các lần chạy trước);
  - **LA** = mọi độ trễ (mặc định của lệnh; `maxlag(1 5)` khi T = 6).
- Các tuỳ chọn khác như trước: `grid(199) refine(4) trim(.15) gridtype(uniform)
  gridsample(effective) history(panel) vce(robust) bwscale(1.5)`; bootstrap
  `boot(499) gridci(100)`, wild. (Nếu 0.4 đổi mặc định, cân nhắc theo.)
- FD/FOD cùng mẫu dữ liệu (cùng seed DGP và dữ liệu thiếu), như trước.
- Coverage chính: `citest(γ0)`; phụ: hull, độ dài, tỷ lệ chạm biên. Không còn
  continuity test.
- Số lần lặp: 400 (Study A), 500 (Study B, kink, power) như trước.

## 2. Danh sách khối

### Study A: ước lượng điểm và SE hệ số (`noboot`)

| Khối | Thiết kế | Ô | Fit | Ghi chú |
|---|---|---:|---:|---|
| A1–A6 | như bản thảo (Bảng 2, Panel A), L3 | 116 | 48.800 | giữ nguyên |
| A7 **mới** | LA; N = 400; T = 6; κ 0, 1; cân bằng/MCAR .30; FD/FOD | 8 | 3.200 | tác động của số độ trễ lên RMSE, SE |
| A-supp | covariance (robust/Windmeijer; jump/kink) như bản thảo | 32 | 12.800 | giữ nguyên |

### Study B: CI của γ và linearity test

| Khối | Thiết kế | Ô | Fit | Ghi chú |
|---|---|---:|---:|---|
| B1, B2, B4, B5 | như bản thảo, L3, bỏ continuity | 40 | 20.000 | coverage + linearity power |
| B3 | linearity size, L3 | 8 | 4.000 | giữ |
| B6 | cấu hình Gong–Seo (FD) | 4 | 2.000 | giữ |
| B7 **mới** | LA; N = 400; κ 0, 1; cân bằng/MCAR .30; FD/FOD | 8 | 4.000 | coverage + linearity ở mặc định |
| B3-LA **mới** | linearity size, LA, N = 400, cân bằng/MCAR .30, FD/FOD | 4 | 2.000 | |
| B8 **mới** | N = 200; L3; κ 0, 1; cân bằng/MCAR .30; FD/FOD | 8 | 4.000 | N nhỏ |
| K | kink ràng buộc (như supplement kink) | 8 | 4.000 | giữ |

### Power của CI (γ0 + c, c = .10/.25/.50; c = 0 lấy từ khối coverage cùng mẫu)

| Khối | Thiết kế | Ô | Fit |
|---|---|---:|---:|
| P-L3 | như Supplement 2 (N 400: cân bằng/gap30/attrition; N 800 cân bằng; κ 0, 1; FD/FOD) | 48 | 24.000 |
| P-LA **mới** | LA; N = 400; cân bằng/gap30; κ 0, 1; FD/FOD | 24 | 12.000 |
| P-GS | cấu hình Gong–Seo, FD (như 2b) | 12 | 6.000 |

**Tổng: 320 ô, khoảng 147.000 fit** (A 64.800; B 36.000; kink 4.000; power 42.000).

## 3. Thời gian ước tính (28 lõi)

Theo thời gian đo được ở các lần chạy trước; khối LA ước chậm hơn L3 khoảng
1,5 lần (nhiều công cụ hơn).

| Phần | Giờ CPU | Giờ máy |
|---|---:|---:|
| Study A (+ A7, A-supp) | ~38 | ~1,5 |
| Study B (B1–B8) | ~160 | ~6 |
| Kink | ~22 | ~1 |
| Power | ~45 | ~2 |
| **Tổng** | **~265** | **~10–12** |

Có thể chạy qua đêm trong một lần; mỗi phần có merge + attestation riêng.

## 4. Bảng dự kiến trong bài

| Bảng | Nội dung | Nguồn |
|---|---|---|
| 3–6 | như bản thảo (Study A) | A1–A6, A-supp |
| mới | L3 và LA: RMSE, coverage hệ số | A7 vs A1 |
| mới | coverage `citest`, hull (độ dài, chạm biên) theo FD/FOD × dữ liệu × κ; hai cấu hình | B1–B8 |
| mới | linearity size và power | B3, B3-LA, B1, B4, B5, B7 |
| mới | power tại γ0 + c, so với Gong–Seo Bảng 2 | P-L3, P-LA, P-GS |
| mới | kink: coverage, độ dài, hệ số | K |

## 5. Gói nộp

`stata/` (lệnh đóng băng) + `mc/` (một harness, registry, log, CSV tổng hợp,
attestation hash) + `sim/` (prototype Python, chỉ dẫn chứng phụ) + README tái
lập (lệnh chạy, thời gian, phiên bản Stata).
