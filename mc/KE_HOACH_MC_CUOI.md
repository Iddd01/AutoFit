# Kế hoạch MC cuối (29/09/2026) — trục chính là thiết kế Gong–Seo

Nguyên tắc: một phiên bản lệnh (0.9.37, đã đóng băng), một harness, một lần
chạy. Ô mốc là thiết kế của Gong–Seo (2026): FD, dữ liệu cân bằng, mọi độ
trễ, lưới quantile 46 điểm p5–p95, trim 0,10, B = 500. Ô mốc phải tái tạo xấp
xỉ Bảng 1–2 của họ (kết quả 2b đã cho thấy điều này). Mọi ô khác chỉ đổi một
yếu tố so với ô mốc: FD → FOD, cân bằng → gap30 / attrition15, N, κ, DGP,
bộ công cụ, T, hoặc mô hình kink. 500 lần lặp mỗi ô; FD và FOD dùng chung mẫu.

Harness: `mc_final_v0937/final/` (xem `FINAL.md`).

## Ước lượng điểm và SE (`final_point_cells.csv`: 198 ô, 99.000 fit)

| Khối | Nội dung | Thay cho |
|---|---|---|
| P1 | N 400/800/1600 × cân bằng/gap30/attr15 × κ 0/.1/.2/.5/1 | Study A khối A1, A4 |
| P2 | N 200 | A2 |
| P3 | q nội sinh, ρ = .9, đuôi t(5) | A3, A5, A6 |
| P4 | `maxlag(1 3)` và `collapse` (T = 6) | mới: tác động của bộ công cụ |
| P5 | T = 10, N 200/400 × ba bộ công cụ | A2/A4 + đánh đổi số công cụ |
| P6 | Windmeijer; mô hình kink (robust/Windmeijer) | supplement của Study A |

## Suy diễn (`final_inf_cells.csv`: 192 ô, 96.000 fit)

| Khối | Nội dung | Thay cho |
|---|---|---|
| I1 | coverage (`citest`, hull) + power linearity; N 400/800 × 3 kiểu dữ liệu × 5 κ | Study B B1, B2, B5 |
| I2 | N 200; q nội sinh | B4 + mới |
| I3 | size của linearity (DGP tuyến tính) | B3 |
| I4 | CI của mô hình kink | supplement kink |
| I5 | power tại γ0 + c (c .10/.25/.50), κ 0/.5/1 | Supplement 2/2b |
| I6 | `maxlag(1 3)`: coverage + power | độ nhạy theo bộ công cụ |

## Thời gian (28 lõi)

Test Stata khoảng 20 phút; POINT khoảng 1–2 giờ; INF khoảng 5 giờ. Tổng
khoảng 7 giờ, chạy qua đêm bằng `mc_final_v0937/run_tonight.ps1`.

## Bảng dự kiến trong bài

| Bảng | Nguồn |
|---|---|
| Tái tạo Gong–Seo (FD cân bằng): coverage và power so với Bảng 1–2 của họ | I1, I5 |
| FD vs FOD × kiểu dữ liệu: γ̂, ρ̂ (RMSE, cặp FD–FOD) | P1, `final_paired.csv` |
| Coverage của hệ số (VCE chung, điều kiện, Windmeijer) | P1, P6 |
| Coverage, độ dài tập tin cậy của γ; power | I1, I5 |
| Linearity: size và power | I3, I1 |
| Độ bền vững: N nhỏ, q nội sinh, ρ = .9, t(5) | P2, P3, I2 |
| Bộ công cụ: mọi độ trễ / 1–3 / collapse; T = 10 | P4, P5, I6 |
| Mô hình kink | P6, I4 |
