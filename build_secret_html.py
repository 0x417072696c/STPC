import csv, numpy as np, matplotlib.pyplot as plt, base64, io

N_list, ln_kf_list = [], []
with open('kineticdb_dataset.csv', 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        try:
            N_list.append(int(row['Length']))
            ln_kf_list.append(float(row['ln_kf']))
        except: pass
N = np.array(N_list); ln_kf = np.array(ln_kf_list)
ln_tau_exp = -ln_kf
tau_exp = np.exp(ln_tau_exp) * 1000

alpha = 0.04
C_opt = np.mean(ln_tau_exp - alpha * N)
ln_tau_stpc_opt = alpha * N + C_opt
tau_stpc_opt = np.exp(ln_tau_stpc_opt) * 1000

coeff = np.polyfit(np.log(N), tau_exp, 1)
tau_log = coeff[0] * np.log(N) + coeff[1]

r_lnkf = np.corrcoef(ln_kf, -np.log(1e-12) - alpha*N)[0,1]
r_lntau = np.corrcoef(ln_tau_exp, ln_tau_stpc_opt)[0,1]
r_tau_stpc = np.corrcoef(np.log(tau_exp), np.log(tau_stpc_opt))[0,1]
r_tau_log  = np.corrcoef(np.log(tau_exp), np.log(np.maximum(tau_log,1e-9)))[0,1]

fig, axes = plt.subplots(3, 1, figsize=(10, 14))
ax = axes[0]
ax.scatter(N, ln_kf, c='#3daee9', s=40)
ax.plot(N, -np.log(1e-12) - alpha*N, '-', color='gold', linewidth=2)
ax.set_xlabel('Длина белка N (а.о.)')
ax.set_ylabel('ln(k_f) (с⁻¹)')
ax.set_title(f'1. Сырая скорость: ln(k_f) vs N (r = {r_lnkf:.4f}, 87 белков)')
ax.grid(alpha=0.3)

ax = axes[1]
ax.scatter(N, ln_tau_exp, c='#3daee9', s=40)
ax.plot(N, ln_tau_stpc_opt, '-', color='gold', linewidth=2, label=f'СТПК (наклон 0.04, смещение {C_opt:.2f})')
ax.set_xlabel('Длина белка N (а.о.)')
ax.set_ylabel('ln(τ) (τ в секундах)')
ax.set_title(f'2. Логарифмическое время: ln(τ) vs N (r = {r_lntau:.4f}, 87 белков)')
ax.legend()
ax.grid(alpha=0.3)

ax = axes[2]
ax.scatter(N, tau_exp, c='#3daee9', s=40, label='KineticDB (87 белков)')
ax.plot(N, tau_stpc_opt, '-', color='gold', linewidth=2, label=f'СТПК (r = {r_tau_stpc:.4f})')
ax.plot(N, tau_log, '--', color='red', linewidth=2, label=f'Логарифм. модель (r = {r_tau_log:.4f})')
ax.set_yscale('log')
ax.set_xlabel('Длина белка N (а.о.)')
ax.set_ylabel('Время сворачивания τ (мс)')
ax.set_title('3. Физическое время: τ (мс) vs N (log-шкала)')
ax.legend()
ax.grid(alpha=0.3)

plt.tight_layout()
buf = io.BytesIO()
plt.savefig(buf, format='png', dpi=150)
buf.seek(0)
img_b64 = base64.b64encode(buf.read()).decode()
buf.close()

html = f'''<!DOCTYPE html>
<html lang="ru">
<head>
<meta charset="UTF-8">
<title>Парадокс Левинталя — решение в модели СТПК</title>
<style>
  body {{ font-family: 'Segoe UI', sans-serif; background: #1a1a1a; color: #eff0f1; padding: 20px; }}
  h1 {{ text-align: center; color: #3daee9; }}
  h2 {{ color: #3daee9; }}
  .formula {{ background: #2d2d2d; padding: 10px; border-radius: 6px; }}
  table {{ border-collapse: collapse; margin: 20px 0; }}
  th, td {{ border: 1px solid #3daee9; padding: 8px; text-align: left; }}
  th {{ background: #333; }}
  img {{ max-width: 100%; margin: 20px 0; }}
  .footer {{ text-align: center; margin-top: 40px; font-size: 0.8rem; color: #555; }}
</style>
</head>
<body>
<h1>Парадокс Левинталя: решение в модели СТПК</h1>
<h2>Фундаментальная формула</h2>
<div class="formula">
  Время сворачивания белка: <b>τ = τ₀ · exp(α·N)</b><br>
  где τ₀ ≈ 10⁻¹² с, α = (Δ/J)² = (0.025 эВ / 0.125 эВ)² = <b>0.04</b>.<br>
  Логарифм времени: <b>ln(τ) = ln(τ₀) + α·N</b> — линейная зависимость.
</div>

<h2>Результаты на независимой выборке KineticDB (87 белков, длины от 36 до 605 а.о.)</h2>
<table>
  <tr><th>Метрика</th><th>Коэффициент корреляции r</th><th>Модель</th></tr>
  <tr><td>ln(k_f) vs N</td><td>{r_lnkf:.4f}</td><td>СТПК</td></tr>
  <tr><td>ln(τ) vs N (оптимальное смещение C={C_opt:.2f})</td><td>{r_lntau:.4f}</td><td>СТПК</td></tr>
  <tr><td>τ (мс) vs N (лог-шкала)</td><td>{r_tau_stpc:.4f}</td><td>СТПК</td></tr>
  <tr><td>τ (мс) vs N (лог-шкала)</td><td>{r_tau_log:.4f}</td><td>Логарифмическая модель</td></tr>
</table>

<img src="data:image/png;base64,{img_b64}" alt="Графики верификации СТПК">

<div class="footer">
  Конфиденциально. © Валентин Лебедкин, 2026. Для служебного пользования.
</div>
</body>
</html>'''

with open('secret_levinthal.html', 'w', encoding='utf-8') as f:
    f.write(html)
print('HTML сохранён: secret_levinthal.html')