'use strict';
/* =====================================================================
   1. LÓGICA MATEMÁTICA  (sin acceso al DOM; verificable por separado)
   ===================================================================== */
// ==== MATH:BEGIN
const SplineMath = (function () {
  const MIN_H = 1e-6; // distancia mínima entre puntos consecutivos (px)

  /** Parámetro t por longitud de cuerda: t0=0, ti = t(i-1) + distancia(Pi-1, Pi) */
  function chordParameters(pts) {
    const t = [0];
    for (let i = 1; i < pts.length; i++) {
      const dx = pts[i].x - pts[i - 1].x, dy = pts[i].y - pts[i - 1].y;
      t.push(t[i - 1] + Math.sqrt(dx * dx + dy * dy));
    }
    return t;
  }

  /** Valida cantidad de puntos y duplicados consecutivos. */
  function validate(pts) {
    const errors = [];
    if (!Array.isArray(pts) || pts.length < 3) {
      errors.push('Se necesitan al menos 3 puntos de control para calcular el spline (recomendado: 8 a 15). Actualmente hay ' + (pts ? pts.length : 0) + '.');
      return { ok: false, errors };
    }
    for (let i = 0; i < pts.length; i++) {
      if (!isFinite(pts[i].x) || !isFinite(pts[i].y)) errors.push('El punto ' + (i + 1) + ' tiene coordenadas no válidas.');
    }
    if (errors.length) return { ok: false, errors };
    for (let i = 1; i < pts.length; i++) {
      const d = Math.hypot(pts[i].x - pts[i - 1].x, pts[i].y - pts[i - 1].y);
      if (d < MIN_H) errors.push('Los puntos consecutivos ' + i + ' y ' + (i + 1) + ' están duplicados (intervalo de parámetro nulo). Muévelos o elimina uno.');
    }
    return { ok: errors.length === 0, errors };
  }

  /** Algoritmo de Thomas para sistema tridiagonal.
      lower[i] multiplica a x[i-1], diag[i] a x[i], upper[i] a x[i+1]. */
  function solveTridiagonal(lower, diag, upper, rhs) {
    const m = diag.length;
    const cp = new Array(m), dp = new Array(m), x = new Array(m);
    let den = diag[0];
    if (Math.abs(den) < 1e-14) throw new Error('Pivote nulo en el sistema tridiagonal.');
    cp[0] = upper[0] / den;
    dp[0] = rhs[0] / den;
    for (let i = 1; i < m; i++) {
      den = diag[i] - lower[i] * cp[i - 1];
      if (Math.abs(den) < 1e-14) throw new Error('Pivote nulo en el sistema tridiagonal.');
      cp[i] = upper[i] / den;
      dp[i] = (rhs[i] - lower[i] * dp[i - 1]) / den;
    }
    x[m - 1] = dp[m - 1];
    for (let i = m - 2; i >= 0; i--) x[i] = dp[i] - cp[i] * x[i + 1];
    return x;
  }

  /** Spline cúbico natural 1D sobre nodos t[0..n] con valores v[0..n]. */
  function naturalSpline1D(t, v) {
    const n = t.length - 1;
    const h = new Array(n);
    for (let i = 0; i < n; i++) h[i] = t[i + 1] - t[i];
    const a = v.slice();
    const c = new Array(n + 1).fill(0); // c0 = cn = 0 (natural)
    if (n >= 2) {
      const m = n - 1; // incógnitas c1..c(n-1)
      const lower = new Array(m).fill(0), diag = new Array(m), upper = new Array(m).fill(0), rhs = new Array(m);
      for (let k = 0; k < m; k++) {
        const i = k + 1;
        lower[k] = h[i - 1];
        diag[k] = 2 * (h[i - 1] + h[i]);
        upper[k] = h[i];
        rhs[k] = 3 * ((a[i + 1] - a[i]) / h[i] - (a[i] - a[i - 1]) / h[i - 1]);
      }
      const sol = solveTridiagonal(lower, diag, upper, rhs);
      for (let k = 0; k < m; k++) c[k + 1] = sol[k];
    }
    const b = new Array(n), d = new Array(n);
    for (let i = 0; i < n; i++) {
      b[i] = (a[i + 1] - a[i]) / h[i] - h[i] * (2 * c[i] + c[i + 1]) / 3;
      d[i] = (c[i + 1] - c[i]) / (3 * h[i]);
    }
    return { a: a.slice(0, n), b, c: c.slice(0, n), d, cEnd: c[n] };
  }

  /** Construye el modelo paramétrico: X(t), Y(t). */
  function build(pts) {
    const v = validate(pts);
    if (!v.ok) throw new Error(v.errors[0]);
    const t = chordParameters(pts);
    const n = pts.length - 1;
    const h = [];
    for (let i = 0; i < n; i++) h.push(t[i + 1] - t[i]);
    const x = naturalSpline1D(t, pts.map(p => p.x));
    const y = naturalSpline1D(t, pts.map(p => p.y));
    return { n, t, h, x, y };
  }

  /** Evalúa el tramo i en tau = t - t_i (forma explícita de la especificación). */
  function evalSegment(model, i, tau) {
    const X = model.x, Y = model.y;
    return {
      x: X.a[i] + X.b[i] * tau + X.c[i] * tau * tau + X.d[i] * tau * tau * tau,
      y: Y.a[i] + Y.b[i] * tau + Y.c[i] * tau * tau + Y.d[i] * tau * tau * tau
    };
  }

  /** Evalúa en un valor global de t. */
  function evalAt(model, t) {
    let i = 0;
    while (i < model.n - 1 && t > model.t[i + 1]) i++;
    return evalSegment(model, i, t - model.t[i]);
  }

  /** Errores numéricos: interpolación, C0, C1, C2 y frontera natural. */
  function verify(model, pts) {
    const n = model.n;
    let eI = 0, e0 = 0, e1 = 0, e2 = 0, eN = 0;
    const comps = [[model.x, 'x'], [model.y, 'y']];
    for (const [co, key] of comps) {
      for (let i = 0; i < n; i++) {
        const h = model.h[i];
        const val = co.a[i] + co.b[i] * h + co.c[i] * h * h + co.d[i] * h * h * h;
        eI = Math.max(eI, Math.abs(co.a[i] - pts[i][key]), Math.abs(val - pts[i + 1][key]));
        if (i < n - 1) {
          e0 = Math.max(e0, Math.abs(val - co.a[i + 1]));
          e1 = Math.max(e1, Math.abs(co.b[i] + 2 * co.c[i] * h + 3 * co.d[i] * h * h - co.b[i + 1]));
          e2 = Math.max(e2, Math.abs(2 * co.c[i] + 6 * co.d[i] * h - 2 * co.c[i + 1]));
        }
      }
      const hl = model.h[n - 1];
      eN = Math.max(eN, Math.abs(2 * co.c[0]), Math.abs(2 * co.c[n - 1] + 6 * co.d[n - 1] * hl));
    }
    return { interpolacion: eI, C0: e0, C1: e1, C2: e2, frontera_natural: eN };
  }

  return { chordParameters, validate, solveTridiagonal, naturalSpline1D, build, evalSegment, evalAt, verify };
})();
// ==== MATH:END

/* =====================================================================
   2. ESTADO GLOBAL Y UTILIDADES
   ===================================================================== */
const $ = id => document.getElementById(id);
const cv = $('cv'), stage = $('stage'), ctx0 = cv.getContext('2d');

const state = {
  img: null, imgName: '', iw: 0, ih: 0,
  points: [],               // coordenadas en píxeles de la imagen ORIGINAL
  model: null, splineOn: false, verify: null,
  mode: 'edit',             // 'edit' | 'anim'
  view: { scale: 1, ox: 0, oy: 0, auto: true },
  cssW: 0, cssH: 0, dpr: 1,
  drag: -1, hover: -1, selected: -1, pan: null, spacePan: false,
  opt: { aux: true, segColors: true, coords: false, fill: false, compare: false, comparePos: 0.5 },
  renderQueued: false, panelsDirty: true
};
const anim = { t: 0, total: 0, phases: [], playing: false, speed: 1, last: 0, raf: 0 };

const SEG_COLORS = ['#ff4d6d', '#ff9f1c', '#2ec4b6', '#4361ee', '#c026d3', '#84cc16'];
const SPLINE_COLOR = '#ff2d6f';
const clamp = (v, a, b) => Math.min(b, Math.max(a, v));
const clamp01 = v => clamp(v, 0, 1);
const easeInOut = p => p < 0.5 ? 2 * p * p : 1 - Math.pow(-2 * p + 2, 2) / 2;
const easeOutBack = p => { const c1 = 1.70158, c3 = c1 + 1; return 1 + c3 * Math.pow(p - 1, 3) + c1 * Math.pow(p - 1, 2); };

function fmt(v, d) {
  if (!isFinite(v)) return '—';
  if (v === 0) return '0';
  const a = Math.abs(v);
  if (a >= 1e5 || a < 1e-3) return v.toExponential(3);
  return v.toFixed(d === undefined ? (a >= 1000 ? 2 : 4) : d);
}
function fmtSig(v) { return isFinite(v) ? String(Number(v.toPrecision(4))) : '—'; }
function mmss(s) { s = Math.max(0, Math.round(s)); return Math.floor(s / 60) + ':' + String(s % 60).padStart(2, '0'); }

function msg(type, text) {
  const el = $('msg');
  el.className = 'msg show ' + type;
  el.textContent = text;
}

/* =====================================================================
   3. VISTA / TRANSFORMACIONES DE COORDENADAS
      Los puntos viven en coordenadas de la imagen original. La vista
      (scale, ox, oy) sólo sirve para mostrarlos en el canvas.
   ===================================================================== */
function fitScale() {
  if (!state.img) return 1;
  const pad = 14;
  return Math.max(0.01, Math.min((state.cssW - 2 * pad) / state.iw, (state.cssH - 2 * pad) / state.ih));
}
function fitView() {
  if (!state.img) { state.view = { scale: 1, ox: 0, oy: 0, auto: true }; return; }
  const s = fitScale();
  state.view.scale = s;
  state.view.ox = (state.cssW - state.iw * s) / 2;
  state.view.oy = (state.cssH - state.ih * s) / 2;
  state.view.auto = true;
}
function zoomAt(factor, cx, cy) {
  if (!state.img) return;
  const v = state.view;
  const ns = clamp(v.scale * factor, fitScale() * 0.3, 64);
  const r = ns / v.scale;
  v.ox = cx - (cx - v.ox) * r;
  v.oy = cy - (cy - v.oy) * r;
  v.scale = ns; v.auto = false;
  scheduleRender();
}
const imgToScreen = (p, v) => ({ x: (v || state.view).ox + p.x * (v || state.view).scale, y: (v || state.view).oy + p.y * (v || state.view).scale });
const screenToImg = (sx, sy) => ({ x: (sx - state.view.ox) / state.view.scale, y: (sy - state.view.oy) / state.view.scale });
function evPos(e) {
  const r = cv.getBoundingClientRect();
  return { x: (e.clientX - r.left) * (state.cssW / r.width), y: (e.clientY - r.top) * (state.cssH / r.height) };
}
function hitTest(sx, sy) {
  let best = -1, bd = 12;
  state.points.forEach((p, i) => {
    const s = imgToScreen(p);
    const d = Math.hypot(s.x - sx, s.y - sy);
    if (d <= bd) { bd = d; best = i; }
  });
  return best;
}

/* =====================================================================
   4. VISUALIZACIÓN (Canvas 2D)
   ===================================================================== */
function roundRectPath(c, x, y, w, h, r) {
  r = Math.min(r, w / 2, h / 2);
  c.beginPath();
  c.moveTo(x + r, y);
  c.arcTo(x + w, y, x + w, y + h, r);
  c.arcTo(x + w, y + h, x, y + h, r);
  c.arcTo(x, y + h, x, y, r);
  c.arcTo(x, y, x + w, y, r);
  c.closePath();
}
function pill(c, text, x, y, ui, bg, fg, anchorRight) {
  c.font = '700 ' + (11 * ui) + 'px system-ui, sans-serif';
  const w = c.measureText(text).width + 9 * ui, h = 16 * ui;
  const px = anchorRight ? x - w : x;
  c.fillStyle = bg; roundRectPath(c, px, y - h, w, h, 4 * ui); c.fill();
  c.fillStyle = fg || '#fff'; c.textBaseline = 'middle'; c.textAlign = 'left';
  c.fillText(text, px + 4.5 * ui, y - h / 2 + 0.5 * ui);
}

/** Puntos (en pantalla) del tramo i entre las fracciones f0..f1 de su intervalo. */
function segPolyline(i, f0, f1, view) {
  const m = state.model, h = m.h[i];
  const N = Math.min(800, Math.max(8, Math.ceil(h * view.scale / 3)));
  const steps = Math.max(2, Math.ceil(N * (f1 - f0)));
  const out = [];
  for (let k = 0; k <= steps; k++) {
    const tau = h * (f0 + (f1 - f0) * k / steps);
    const p = SplineMath.evalSegment(m, i, tau);
    out.push({ x: view.ox + p.x * view.scale, y: view.oy + p.y * view.scale });
  }
  return out;
}
function strokePoly(c, pts, color, width, dash) {
  if (pts.length < 2) return;
  c.beginPath(); c.moveTo(pts[0].x, pts[0].y);
  for (let k = 1; k < pts.length; k++) c.lineTo(pts[k].x, pts[k].y);
  c.strokeStyle = color; c.lineWidth = width; c.setLineDash(dash || []);
  c.lineJoin = 'round'; c.lineCap = 'round'; c.stroke(); c.setLineDash([]);
}
function zoneRect() {
  const P = state.points;
  let x0 = Infinity, y0 = Infinity, x1 = -Infinity, y1 = -Infinity;
  P.forEach(p => { x0 = Math.min(x0, p.x); y0 = Math.min(y0, p.y); x1 = Math.max(x1, p.x); y1 = Math.max(y1, p.y); });
  const pad = 0.12 * Math.max(x1 - x0, y1 - y0, 20) + 10;
  return { x0: Math.max(0, x0 - pad), y0: Math.max(0, y0 - pad), x1: Math.min(state.iw, x1 + pad), y1: Math.min(state.ih, y1 + pad) };
}

/** Dibuja toda la escena. Se usa en pantalla (ui=1) y en la exportación PNG (ui>=1). */
function drawScene(c, W, H, view, ui, o) {
  c.save();
  c.fillStyle = '#10141b'; c.fillRect(0, 0, W, H);
  if (!state.img) { c.restore(); return; }
  const iw = state.iw, ih = state.ih, S = view.scale;
  c.imageSmoothingEnabled = true; c.imageSmoothingQuality = 'high';
  c.drawImage(state.img, view.ox, view.oy, iw * S, ih * S);

  const model = state.model;
  const P = state.points;
  const sp = P.map(p => imgToScreen(p, view));
  const n = P.length;

  // --- Capas que quedan a la derecha del divisor de comparación ---
  const cmp = (o.compare !== null && o.compare !== undefined) ? o.compare : null;
  const divX = cmp !== null ? view.ox + cmp * iw * S : 0;
  c.save();
  if (cmp !== null) { c.beginPath(); c.rect(divX, 0, Math.max(0, W - divX), H); c.clip(); }

  // Zona seleccionada (oscurece el exterior)
  if (o.zone > 0 && n >= 1) {
    const z = zoneRect();
    const zx = view.ox + z.x0 * S, zy = view.oy + z.y0 * S, zw = (z.x1 - z.x0) * S, zh = (z.y1 - z.y0) * S;
    c.save();
    c.beginPath();
    c.rect(view.ox, view.oy, iw * S, ih * S);
    c.rect(zx, zy, zw, zh);
    c.fillStyle = 'rgba(8,12,24,' + (0.58 * o.zone) + ')';
    c.fill('evenodd');
    c.globalAlpha = o.zone;
    c.setLineDash([10 * ui, 6 * ui]); c.lineWidth = 2.4 * ui; c.strokeStyle = '#fde047';
    c.strokeRect(zx, zy, zw, zh); c.setLineDash([]);
    pill(c, 'ZONA SELECCIONADA', zx, zy - 4 * ui, ui, '#a16207', '#fff');
    c.restore();
  }

  // Relleno del perfil reconstruido (polígono cerrado por la cuerda final)
  if (o.fill > 0 && model) {
    let all = [];
    for (let i = 0; i < model.n; i++) { const pl = segPolyline(i, 0, 1, view); all = all.concat(i ? pl.slice(1) : pl); }
    c.save();
    c.beginPath(); c.moveTo(all[0].x, all[0].y);
    for (let k = 1; k < all.length; k++) c.lineTo(all[k].x, all[k].y);
    c.closePath();
    c.fillStyle = 'rgba(255,170,40,' + (0.42 * o.fill) + ')'; c.fill();
    c.globalAlpha = 0.85 * o.fill;
    strokePoly(c, [all[all.length - 1], all[0]], '#ffffff', 2 * ui, [7 * ui, 6 * ui]);
    c.restore();
  }

  // Líneas auxiliares
  const shown = Math.min(n, Math.ceil(o.points - 1e-9));
  const fullCnt = Math.floor(o.points + 1e-9);
  if (o.aux && fullCnt >= 2) {
    c.save(); c.globalAlpha = (o.auxAlpha === undefined ? 0.9 : o.auxAlpha);
    const pl = sp.slice(0, Math.min(fullCnt, n));
    strokePoly(c, pl, 'rgba(0,0,0,0.55)', 3.6 * ui);
    strokePoly(c, pl, '#ffe066', 1.7 * ui, [6 * ui, 5 * ui]);
    c.restore();
  }

  // Segmentos del spline
  let tip = null;
  if (model && o.segs > 0) {
    const polys = [];
    for (let i = 0; i < model.n; i++) {
      const f = Math.min(1, o.segs - i);
      if (f <= 1e-9) break;
      const pl = segPolyline(i, 0, f, view);
      polys.push({ i, pl });
      tip = pl[pl.length - 1];
    }
    polys.forEach(s => strokePoly(c, s.pl, 'rgba(0,0,0,0.6)', (s.i === o.highlight ? 7.5 : 6) * ui));
    polys.forEach(s => strokePoly(c, s.pl, o.segColors ? SEG_COLORS[s.i % SEG_COLORS.length] : SPLINE_COLOR, (s.i === o.highlight ? 4.6 : 3.4) * ui));
    if (o.tip && o.segs < model.n - 1e-9 && tip) {
      c.beginPath(); c.arc(tip.x, tip.y, 8 * ui, 0, Math.PI * 2);
      c.fillStyle = 'rgba(253,224,71,0.35)'; c.fill();
      c.beginPath(); c.arc(tip.x, tip.y, 4 * ui, 0, Math.PI * 2);
      c.fillStyle = '#fde047'; c.fill(); c.lineWidth = 1.5 * ui; c.strokeStyle = '#000'; c.stroke();
    }
  }

  // Puntos de control con numeración
  for (let i = 0; i < shown; i++) {
    const frac = clamp01(o.points - i);
    const sc = frac < 1 ? Math.max(0.01, easeOutBack(frac)) : 1;
    const isSel = (i === o.selected), isHov = (i === o.hover);
    const r = (isSel ? 9 : isHov ? 8.5 : 7) * ui * sc;
    const p = sp[i];
    c.beginPath(); c.arc(p.x, p.y, r + 1.6 * ui, 0, Math.PI * 2); c.strokeStyle = 'rgba(0,0,0,0.75)'; c.lineWidth = 3 * ui; c.stroke();
    c.beginPath(); c.arc(p.x, p.y, r, 0, Math.PI * 2);
    c.strokeStyle = isSel ? '#fde047' : isHov ? '#67e8f9' : '#ffffff'; c.lineWidth = 2.2 * ui; c.stroke();
    c.beginPath(); c.arc(p.x, p.y, 2.3 * ui * sc, 0, Math.PI * 2); c.fillStyle = '#2563eb'; c.fill();
    c.lineWidth = 1 * ui; c.strokeStyle = '#fff'; c.stroke();
    const bg = (i === 0 && n > 1) ? '#15803d' : (i === n - 1 && n > 1) ? '#b91c1c' : 'rgba(15,23,42,0.9)';
    let label = String(i + 1);
    if (o.coords) label += '  (' + Math.round(P[i].x) + ', ' + Math.round(P[i].y) + ')';
    c.save(); c.globalAlpha = clamp01(frac * 1.5);
    pill(c, label, p.x + 9 * ui, p.y - 8 * ui, ui, bg, '#fff');
    c.restore();
  }
  c.restore(); // fin del clip de comparación

  // Divisor de comparación antes / después
  if (cmp !== null && cmp > 0.002) {
    c.save();
    c.beginPath(); c.moveTo(divX, 0); c.lineTo(divX, H);
    c.strokeStyle = 'rgba(0,0,0,0.5)'; c.lineWidth = 5 * ui; c.stroke();
    c.strokeStyle = '#fff'; c.lineWidth = 2.2 * ui; c.stroke();
    c.beginPath(); c.arc(divX, H / 2, 13 * ui, 0, Math.PI * 2); c.fillStyle = '#fff'; c.fill();
    c.fillStyle = '#0f172a'; c.font = '700 ' + (13 * ui) + 'px system-ui'; c.textAlign = 'center'; c.textBaseline = 'middle';
    c.fillText('↔', divX, H / 2 + 0.5 * ui);
    pill(c, 'ORIGINAL', divX - 8 * ui, 30 * ui, ui, 'rgba(15,23,42,0.9)', '#fff', true);
    pill(c, 'RECONSTRUCCIÓN', divX + 8 * ui, 30 * ui, ui, '#c2410c', '#fff', false);
    c.restore();
  }
  c.restore();
}

/** Rótulos de la animación (se dibujan sólo en pantalla). */
function drawCaption(c, W, H, cap) {
  if (!cap) return;
  c.save();
  c.font = '700 13px system-ui, sans-serif'; c.textBaseline = 'middle'; c.textAlign = 'left';
  const head = cap.step + '  ·  ' + cap.title.toUpperCase();
  const hw = c.measureText(head).width + 28;
  c.fillStyle = 'rgba(15,23,42,0.86)'; roundRectPath(c, 14, 14, hw, 30, 8); c.fill();
  c.fillStyle = '#fff'; c.fillText(head, 28, 29.5);

  const lines = [{ t: cap.desc, f: '600 14px system-ui, sans-serif', col: '#f8fafc' }]
    .concat((cap.eq || []).map(t => ({ t, f: '12px ui-monospace, Menlo, Consolas, monospace', col: '#fde68a' })));
  let bw = 0;
  lines.forEach(l => { c.font = l.f; bw = Math.max(bw, c.measureText(l.t).width); });
  bw = Math.min(W - 28, bw + 28);
  const bh = 14 + lines.length * 21;
  c.fillStyle = 'rgba(15,23,42,0.86)'; roundRectPath(c, 14, H - 14 - bh, bw, bh, 10); c.fill();
  lines.forEach((l, k) => { c.font = l.f; c.fillStyle = l.col; c.fillText(l.t, 28, H - 14 - bh + 17 + k * 21, bw - 24); });
  c.restore();
}

/* ---------- Opciones de dibujo ---------- */
function editOptions() {
  const o = state.opt;
  return {
    points: state.points.length, aux: o.aux, auxAlpha: 0.9,
    segs: state.model ? state.model.n : 0, fill: (o.fill && state.model) ? 1 : 0, zone: 0,
    segColors: o.segColors, coords: o.coords, highlight: -1, tip: false,
    hover: state.hover, selected: state.selected,
    compare: o.compare ? o.comparePos : null, caption: null
  };
}
function exportOptions() {
  const o = editOptions(); o.hover = -1; o.selected = -1; o.compare = null; return o;
}

/* ---------- Render principal ---------- */
function scheduleRender(panels) {
  if (panels) state.panelsDirty = true;
  if (state.renderQueued) return;
  state.renderQueued = true;
  requestAnimationFrame(() => {
    state.renderQueued = false;
    render();
    if (state.panelsDirty) { state.panelsDirty = false; updatePanels(); }
  });
}
function render() {
  const c = ctx0, W = state.cssW, H = state.cssH;
  c.setTransform(state.dpr, 0, 0, state.dpr, 0, 0);
  const o = state.mode === 'anim' ? animOptions() : editOptions();
  drawScene(c, W, H, state.view, 1, o);
  if (state.mode === 'anim') drawCaption(c, W, H, o.caption);
  $('zoomLbl').textContent = state.img ? Math.round(state.view.scale * 100) + '%' : '—';
  $('emptyHint').style.display = state.img ? 'none' : 'flex';
}
function resizeCanvas() {
  const r = stage.getBoundingClientRect();
  const w = Math.max(100, Math.round(r.width)), h = Math.max(100, Math.round(r.height));
  state.cssW = w; state.cssH = h; state.dpr = window.devicePixelRatio || 1;
  cv.width = Math.round(w * state.dpr); cv.height = Math.round(h * state.dpr);
  cv.style.width = w + 'px'; cv.style.height = h + 'px';
  if (state.view.auto) fitView();
  scheduleRender();
}

/* =====================================================================
   5. GESTIÓN DE PUNTOS Y MODELO
   ===================================================================== */
function refreshModel() {
  if (!state.splineOn) { state.model = null; state.verify = null; return; }
  const v = SplineMath.validate(state.points);
  if (!v.ok) {
    state.model = null; state.verify = null;
    msg('warn', 'El spline se desactivó temporalmente: ' + v.errors[0]);
    return;
  }
  try {
    state.model = SplineMath.build(state.points);
    state.verify = SplineMath.verify(state.model, state.points);
  } catch (e) {
    state.model = null; state.verify = null;
    msg('error', 'Error numérico al calcular el spline: ' + e.message);
  }
}
function afterPointsChanged() {
  if (state.points.length === 0) state.splineOn = false;
  refreshModel();
  scheduleRender(true);
}
function addPoint(p) {
  const last = state.points[state.points.length - 1];
  if (last && Math.hypot(p.x - last.x, p.y - last.y) < 0.5) {
    msg('warn', 'Ese punto coincide con el anterior (produciría un intervalo de parámetro nulo). Elige otra posición.');
    return false;
  }
  state.points.push({ x: p.x, y: p.y });
  state.selected = state.points.length - 1;
  afterPointsChanged();
  return true;
}
function removeLast() {
  if (state.mode !== 'edit') return;
  if (!state.points.length) { msg('info', 'No hay puntos para eliminar.'); return; }
  state.points.pop(); state.selected = -1; state.hover = -1;
  afterPointsChanged();
}
function clearPoints() {
  if (state.mode !== 'edit') return;
  state.points = []; state.selected = -1; state.hover = -1;
  state.splineOn = false; state.model = null; state.verify = null;
  msg('info', 'Se eliminaron todos los puntos.');
  afterPointsChanged();
}
function calculateSpline() {
  if (!state.img) { msg('error', 'Primero carga una fotografía y coloca los puntos de control sobre ella.'); return false; }
  const v = SplineMath.validate(state.points);
  if (!v.ok) { msg('error', v.errors.join(' ')); return false; }
  try {
    state.model = SplineMath.build(state.points);
    state.verify = SplineMath.verify(state.model, state.points);
    state.splineOn = true;
  } catch (e) {
    state.model = null; state.splineOn = false;
    msg('error', 'Error numérico al calcular el spline: ' + e.message);
    scheduleRender(true);
    return false;
  }
  const n = state.points.length;
  let text = 'Spline calculado: ' + state.model.n + ' segmentos cúbicos para ' + n + ' puntos. La curva pasa exactamente por todos los puntos de control.';
  let type = 'ok';
  if (n < 8 || n > 15) { text += ' Nota: se recomiendan entre 8 y 15 puntos.'; type = 'info'; }
  msg(type, text);
  scheduleRender(true);
  return true;
}

/* =====================================================================
   6. INTERACCIÓN (ratón / teclado / archivos)
   ===================================================================== */
function updateReadout(sx, sy) {
  if (!state.img) { $('coordReadout').textContent = 'Cursor: —'; return; }
  const p = screenToImg(sx, sy);
  const inside = p.x >= 0 && p.y >= 0 && p.x <= state.iw && p.y <= state.ih;
  $('coordReadout').textContent = inside
    ? 'Cursor: x = ' + p.x.toFixed(1) + ' · y = ' + p.y.toFixed(1) + ' px (imagen ' + state.iw + '×' + state.ih + ')'
    : 'Cursor: fuera de la imagen';
}

cv.addEventListener('mousedown', e => { if (e.button === 1) e.preventDefault(); });
cv.addEventListener('contextmenu', e => {
  e.preventDefault();
  if (state.mode !== 'edit' || !state.img) return;
  const { x, y } = evPos(e);
  const h = hitTest(x, y);
  if (h >= 0) {
    state.points.splice(h, 1); state.selected = -1; state.hover = -1;
    afterPointsChanged();
  }
});
cv.addEventListener('pointerdown', e => {
  if (state.mode !== 'edit') return;
  const { x: sx, y: sy } = evPos(e);
  if (e.button === 1 || (e.button === 0 && (e.shiftKey || state.spacePan))) {
    state.pan = { sx, sy, ox: state.view.ox, oy: state.view.oy };
    cv.setPointerCapture(e.pointerId); cv.style.cursor = 'grabbing'; e.preventDefault(); return;
  }
  if (e.button !== 0) return;
  if (!state.img) { msg('error', 'Primero carga una fotografía con el botón “Cargar fotografía”.'); return; }
  const hit = hitTest(sx, sy);
  if (hit >= 0) {
    state.drag = hit; state.selected = hit;
    cv.setPointerCapture(e.pointerId); scheduleRender(true); return;
  }
  const p = screenToImg(sx, sy);
  if (p.x < 0 || p.y < 0 || p.x > state.iw || p.y > state.ih) { msg('warn', 'Haz clic dentro de la fotografía.'); return; }
  addPoint(p);
});
cv.addEventListener('pointermove', e => {
  const { x: sx, y: sy } = evPos(e);
  updateReadout(sx, sy);
  if (state.pan) {
    state.view.ox = state.pan.ox + (sx - state.pan.sx);
    state.view.oy = state.pan.oy + (sy - state.pan.sy);
    state.view.auto = false; scheduleRender(); return;
  }
  if (state.drag >= 0) {
    const p = screenToImg(sx, sy);
    state.points[state.drag] = { x: clamp(p.x, 0, state.iw), y: clamp(p.y, 0, state.ih) };
    afterPointsChanged(); return;
  }
  if (state.mode === 'edit') {
    const h = hitTest(sx, sy);
    cv.style.cursor = h >= 0 ? 'grab' : (state.img ? 'crosshair' : 'default');
    if (h !== state.hover) { state.hover = h; scheduleRender(); }
  }
});
function endPointer(e) {
  if (state.pan) { state.pan = null; cv.style.cursor = state.img ? 'crosshair' : 'default'; }
  if (state.drag >= 0) { state.drag = -1; scheduleRender(true); }
  try { cv.releasePointerCapture(e.pointerId); } catch (_) {}
}
cv.addEventListener('pointerup', endPointer);
cv.addEventListener('pointercancel', endPointer);
cv.addEventListener('pointerleave', () => {
  $('coordReadout').textContent = 'Cursor: —';
  if (state.hover !== -1 && state.drag < 0) { state.hover = -1; scheduleRender(); }
});
cv.addEventListener('wheel', e => {
  if (!state.img) return;
  e.preventDefault();
  const { x, y } = evPos(e);
  zoomAt(e.deltaY < 0 ? 1.15 : 1 / 1.15, x, y);
}, { passive: false });

window.addEventListener('keydown', e => {
  const tag = (e.target.tagName || '').toLowerCase();
  if (['input', 'textarea', 'select', 'button'].includes(tag)) return;
  if (e.code === 'Space') { state.spacePan = true; e.preventDefault(); return; }
  if (state.mode !== 'edit') return;
  if (e.key === 'Delete' || e.key === 'Backspace' || (e.ctrlKey && e.key.toLowerCase() === 'z')) {
    removeLast(); e.preventDefault(); return;
  }
  const s = state.selected;
  if (s >= 0 && state.points[s] && e.key.startsWith('Arrow')) {
    const st = e.shiftKey ? 5 : 1;
    const p = state.points[s];
    if (e.key === 'ArrowLeft') p.x -= st;
    if (e.key === 'ArrowRight') p.x += st;
    if (e.key === 'ArrowUp') p.y -= st;
    if (e.key === 'ArrowDown') p.y += st;
    p.x = clamp(p.x, 0, state.iw); p.y = clamp(p.y, 0, state.ih);
    afterPointsChanged(); e.preventDefault();
  }
});
window.addEventListener('keyup', e => { if (e.code === 'Space') state.spacePan = false; });

function loadImageFile(file) {
  if (!file) return;
  if (!file.type || !file.type.startsWith('image/')) { msg('error', 'El archivo seleccionado no es una imagen válida.'); return; }
  const url = URL.createObjectURL(file);
  const im = new Image();
  im.onload = () => {
    state.img = im; state.imgName = file.name; state.iw = im.naturalWidth; state.ih = im.naturalHeight;
    state.points = []; state.model = null; state.verify = null; state.splineOn = false;
    state.selected = -1; state.hover = -1;
    fitView(); URL.revokeObjectURL(url);
    msg('info', 'Imagen cargada: ' + file.name + ' (' + state.iw + '×' + state.ih + ' px). Haz clic sobre el contorno para colocar los puntos.');
    cv.style.cursor = 'crosshair';
    scheduleRender(true);
  };
  im.onerror = () => { URL.revokeObjectURL(url); msg('error', 'No se pudo leer la imagen. Prueba con un archivo JPG, PNG o WebP.'); };
  im.src = url;
}
$('fileInput').addEventListener('change', e => { loadImageFile(e.target.files[0]); e.target.value = ''; });
stage.addEventListener('dragover', e => { e.preventDefault(); });
stage.addEventListener('drop', e => {
  e.preventDefault();
  if (state.mode !== 'edit') return;
  const f = e.dataTransfer.files && e.dataTransfer.files[0];
  if (f) loadImageFile(f);
});

$('undoBtn').addEventListener('click', removeLast);
$('clearBtn').addEventListener('click', clearPoints);
$('calcBtn').addEventListener('click', () => { if (state.mode === 'anim') setMode('edit'); calculateSpline(); });
$('zoomInBtn').addEventListener('click', () => zoomAt(1.3, state.cssW / 2, state.cssH / 2));
$('zoomOutBtn').addEventListener('click', () => zoomAt(1 / 1.3, state.cssW / 2, state.cssH / 2));
$('zoomFitBtn').addEventListener('click', () => { fitView(); scheduleRender(); });
$('optAux').addEventListener('change', e => { state.opt.aux = e.target.checked; scheduleRender(); });
$('optSeg').addEventListener('change', e => { state.opt.segColors = e.target.checked; scheduleRender(); });
$('optCoords').addEventListener('change', e => { state.opt.coords = e.target.checked; scheduleRender(); });
$('optFill').addEventListener('change', e => { state.opt.fill = e.target.checked; scheduleRender(); });
$('optCompare').addEventListener('change', e => { state.opt.compare = e.target.checked; scheduleRender(); });
$('comparePos').addEventListener('input', e => { state.opt.comparePos = e.target.value / 100; if (state.opt.compare) scheduleRender(); });

/* =====================================================================
   7. ANIMACIÓN (usa los coeficientes reales del spline calculado)
      Pasos: 1 Pieza original · 2 Zona seleccionada · 3 Puntos de control
             4 Primeros segmentos · 5 Unión (spline completo)
             6 Reconstrucción digital · 7 Resultado final
   ===================================================================== */
const PHASE_NAMES = ['Pieza original', 'Zona seleccionada', 'Puntos de control', 'Primeros segmentos',
  'Unión de segmentos (spline completo)', 'Reconstrucción digital', 'Resultado final'];

function buildPhases() {
  const n = state.points.length, ns = state.model.n;
  const durs = [2.0, 2.5, Math.min(6, Math.max(3, n * 0.5)), 3.0,
    Math.min(10, Math.max(3, (ns - 1) * 1.0 + 1.5)), 3.5, 6.0];
  let s = 0;
  anim.phases = PHASE_NAMES.map((name, k) => { const p = { name, dur: durs[k], start: s }; s += durs[k]; return p; });
  anim.total = s;
  $('steps').innerHTML = anim.phases.map((p, k) => '<button class="chip" data-k="' + k + '">' + (k + 1) + '. ' + p.name + '</button>').join('');
}
function segEquations(i) {
  const m = state.model, ti = m.t[i];
  const poly = (nm, co) => nm + i + '(t) = ' + fmtSig(co.a[i]) + ' + ' + fmtSig(co.b[i]) + '·τ + ' + fmtSig(co.c[i]) + '·τ² + ' + fmtSig(co.d[i]) + '·τ³';
  return [
    'Tramo ' + (i + 1) + ' de ' + m.n + ':  t ∈ [' + fmtSig(ti) + ', ' + fmtSig(m.t[i + 1]) + ']   τ = t − ' + fmtSig(ti),
    poly('X', m.x), poly('Y', m.y)
  ];
}
function animOptions() {
  const m = state.model, ns = m ? m.n : 0, n = state.points.length;
  const t = anim.t;
  let k = anim.phases.length - 1;
  for (let i = 0; i < anim.phases.length; i++) if (t < anim.phases[i].start + anim.phases[i].dur) { k = i; break; }
  const ph = anim.phases[k], p = clamp01((t - ph.start) / ph.dur);
  const o = {
    points: 0, aux: false, auxAlpha: 0.8, segs: 0, fill: 0, zone: 0, compare: null,
    segColors: state.opt.segColors, coords: false, highlight: -1, tip: true, hover: -1, selected: -1, phase: k, p
  };
  const cap = { step: 'PASO ' + (k + 1) + '/7', title: ph.name, desc: '', eq: [] };
  switch (k) {
    case 0:
      cap.desc = 'Fotografía original de la pieza arqueológica'; break;
    case 1:
      o.zone = easeInOut(p); cap.desc = 'Zona del contorno que se va a reconstruir'; break;
    case 2:
      o.zone = 1; o.points = p * n; o.aux = true; o.coords = true;
      cap.desc = 'Puntos de control Pᵢ = (xᵢ, yᵢ): ' + Math.min(n, Math.ceil(o.points - 1e-9)) + ' de ' + n + ' colocados sobre el perfil'; break;
    case 3:
      o.zone = 1; o.points = n; o.aux = true; o.coords = true; o.segs = easeInOut(p); o.highlight = 0;
      cap.desc = 'Primer tramo: polinomio cúbico paramétrico X₀(t), Y₀(t)'; cap.eq = segEquations(0); break;
    case 4: {
      o.zone = 1; o.points = n; o.aux = true; o.segs = 1 + (ns - 1) * easeInOut(Math.min(1, p / 0.8));
      const cur = clamp(Math.ceil(o.segs - 1e-9) - 1, 0, ns - 1);
      if (p < 0.8) { cap.desc = 'Se incorporan los tramos restantes (segmento ' + (cur + 1) + ' de ' + ns + ')'; cap.eq = segEquations(cur); }
      else { cap.desc = 'Spline completo: ' + ns + ' tramos unidos con continuidad C⁰, C¹ y C²'; }
      break;
    }
    case 5:
      o.zone = 1 - 0.6 * easeInOut(p); o.points = n; o.aux = true; o.auxAlpha = 0.5 * (1 - p);
      o.segs = ns; o.fill = easeInOut(p);
      cap.desc = 'Relleno del perfil generado a partir del spline calculado'; break;
    case 6: {
      o.points = n; o.segs = ns; o.fill = 1; o.zone = 0.4 * (1 - easeInOut(Math.min(1, p / 0.2)));
      let cmp = 0;
      if (p > 0.12) { const q = clamp01((p - 0.12) / 0.76); cmp = q < 0.5 ? easeInOut(q * 2) : easeInOut(2 - q * 2); }
      o.compare = cmp;
      cap.desc = p < 0.9 ? 'Comparación antes / después sobre la fotografía original' : 'Resultado final: perfil reconstruido mediante spline cúbico paramétrico';
      break;
    }
  }
  o.caption = cap;
  return o;
}
function setMode(m) {
  state.mode = m;
  document.body.classList.toggle('anim-mode', m === 'anim');
  document.querySelectorAll('.edit-only').forEach(el => {
    if ('disabled' in el) el.disabled = (m === 'anim');
    el.classList.toggle('disabled', m === 'anim');
  });
  if (m === 'edit') { anim.playing = false; cancelAnimationFrame(anim.raf); }
  $('scrub').disabled = (m !== 'anim');
  $('exitAnimBtn').disabled = (m !== 'anim');
  $('resetBtn').disabled = (m !== 'anim');
  updateAnimUI();
  scheduleRender();
}
function startAnimation() {
  if (!state.img) { msg('error', 'Primero carga una fotografía y coloca los puntos de control.'); return; }
  if (!state.model && !calculateSpline()) return;
  if (state.mode === 'edit') { buildPhases(); fitView(); setMode('anim'); anim.t = 0; }
  playAnim();
}
function playAnim() {
  if (state.mode !== 'anim') { startAnimation(); return; }
  if (anim.t >= anim.total - 1e-6) anim.t = 0;
  anim.playing = true; anim.last = performance.now();
  cancelAnimationFrame(anim.raf); anim.raf = requestAnimationFrame(tick);
  updateAnimUI();
}
function pauseAnim() { anim.playing = false; cancelAnimationFrame(anim.raf); updateAnimUI(); }
function resetAnim() {
  if (state.mode !== 'anim') return;
  anim.playing = false; cancelAnimationFrame(anim.raf); anim.t = 0;
  updateAnimUI(); scheduleRender();
}
function tick(ts) {
  if (!anim.playing) return;
  const dt = Math.min(0.1, (ts - anim.last) / 1000);
  anim.last = ts;
  anim.t = Math.min(anim.total, anim.t + dt * anim.speed);
  render(); updateAnimUI();
  if (anim.t >= anim.total) { anim.playing = false; updateAnimUI(); return; }
  anim.raf = requestAnimationFrame(tick);
}
function updateAnimUI() {
  const inAnim = state.mode === 'anim';
  $('playBtn').disabled = anim.playing;
  $('pauseBtn').disabled = !anim.playing;
  $('scrub').value = inAnim && anim.total ? Math.round(anim.t / anim.total * 1000) : 0;
  $('timeLbl').textContent = inAnim ? mmss(anim.t) + ' / ' + mmss(anim.total) : '0:00 / 0:00';
  let k = -1;
  if (inAnim) { k = anim.phases.length - 1; for (let i = 0; i < anim.phases.length; i++) if (anim.t < anim.phases[i].start + anim.phases[i].dur) { k = i; break; } }
  document.querySelectorAll('#steps .chip').forEach((el, i) => {
    el.classList.toggle('active', i === k); el.classList.toggle('done', k >= 0 && i < k); el.disabled = !inAnim;
  });
}
$('animBtn').addEventListener('click', () => { if (state.mode === 'anim') { anim.t = 0; playAnim(); } else startAnimation(); });
$('playBtn').addEventListener('click', playAnim);
$('pauseBtn').addEventListener('click', pauseAnim);
$('resetBtn').addEventListener('click', resetAnim);
$('exitAnimBtn').addEventListener('click', () => { setMode('edit'); msg('info', 'Modo edición: puedes mover los puntos y recalcular el spline.'); });
$('speedRange').addEventListener('input', e => { anim.speed = parseFloat(e.target.value); $('speedLbl').textContent = anim.speed.toFixed(2) + '×'; });
$('scrub').addEventListener('input', e => { anim.t = e.target.value / 1000 * anim.total; updateAnimUI(); scheduleRender(); });
$('steps').addEventListener('click', e => {
  const b = e.target.closest('.chip'); if (!b || state.mode !== 'anim') return;
  const k = parseInt(b.dataset.k, 10);
  anim.playing = false; cancelAnimationFrame(anim.raf);
  // Salta al final del paso elegido (muestra ese paso ya completo)
  anim.t = k === 0 ? 0 : anim.phases[k].start + anim.phases[k].dur * (k === 6 ? 0.97 : 0.999);
  updateAnimUI(); scheduleRender();
});

/* =====================================================================
   8. EXPORTACIÓN (PNG y JSON)
   ===================================================================== */
function downloadBlob(blob, name) {
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url; a.download = name; document.body.appendChild(a); a.click();
  setTimeout(() => { URL.revokeObjectURL(url); a.remove(); }, 500);
}
function baseName() { return (state.imgName || 'pieza').replace(/\.[^.]+$/, ''); }

function exportPNG() {
  if (!state.img) { msg('error', 'No hay fotografía cargada para exportar.'); return; }
  const L = Math.max(state.iw, state.ih);
  const k = Math.min(4096 / L, Math.max(1, 1920 / L));
  const W = Math.round(state.iw * k), H = Math.round(state.ih * k);
  const off = document.createElement('canvas'); off.width = W; off.height = H;
  const c = off.getContext('2d');
  const ui = Math.max(1, Math.max(W, H) / 900);
  drawScene(c, W, H, { scale: k, ox: 0, oy: 0 }, ui, exportOptions());
  off.toBlob(b => {
    if (!b) { msg('error', 'No se pudo generar la imagen PNG.'); return; }
    downloadBlob(b, baseName() + '_spline.png');
    msg('ok', 'PNG exportado (' + W + '×' + H + ' px)' + (state.model ? ' con puntos numerados y spline.' : ' con los puntos. Calcula el spline para incluir la curva.'));
  }, 'image/png');
}
function exportJSON() {
  if (!state.model) { msg('error', 'Primero calcula el spline (botón CALCULAR SPLINE) para exportar los resultados.'); return; }
  const m = state.model, P = state.points;
  const data = {
    aplicacion: 'Reconstrucción arqueológica mediante Spline Cúbico Paramétrico Natural',
    fecha: new Date().toISOString(),
    imagen: {
      nombre: state.imgName, ancho_px: state.iw, alto_px: state.ih,
      sistema_coordenadas: 'Píxeles de la imagen original; origen en la esquina superior izquierda; Y crece hacia abajo'
    },
    configuracion_modelo: {
      metodo: 'Spline cúbico paramétrico natural',
      parametrizacion: 'Longitud de cuerda: t0 = 0, ti = t(i-1) + sqrt((xi-x(i-1))^2 + (yi-y(i-1))^2)',
      polinomio: 'X_i(t) = a + b(t-t_i) + c(t-t_i)^2 + d(t-t_i)^3 ; igual para Y_i(t); t en [t_i, t_(i+1)]',
      condicion_frontera: 'Natural: X\'\'(t0) = X\'\'(tn) = 0 y Y\'\'(t0) = Y\'\'(tn) = 0',
      solucion_sistema: 'Sistema tridiagonal resuelto con el algoritmo de Thomas (X e Y por separado)',
      cantidad_puntos: P.length,
      cantidad_intervalos: m.n,
      nota_indices: 'puntos[].indice empieza en 1; intervalos[].i empieza en 0 (como en las fórmulas)'
    },
    puntos: P.map((p, i) => ({ indice: i + 1, x: p.x, y: p.y, t: m.t[i] })),
    parametros_t: m.t.slice(),
    intervalos: m.h.map((h, i) => ({
      i, t_inicio: m.t[i], t_fin: m.t[i + 1], h,
      X: { a: m.x.a[i], b: m.x.b[i], c: m.x.c[i], d: m.x.d[i] },
      Y: { a: m.y.a[i], b: m.y.b[i], c: m.y.c[i], d: m.y.d[i] }
    })),
    verificacion_numerica_px: state.verify
  };
  downloadBlob(new Blob([JSON.stringify(data, null, 2)], { type: 'application/json' }), baseName() + '_spline.json');
  msg('ok', 'JSON exportado con coordenadas, parámetros t, coeficientes y configuración del modelo.');
}
$('pngBtn').addEventListener('click', exportPNG);
$('jsonBtn').addEventListener('click', exportJSON);

/* =====================================================================
   9. PANELES INFORMATIVOS (tablas y resumen)
   ===================================================================== */
function updatePanels() {
  const P = state.points, m = state.model;
  const n = P.length;
  $('countLbl').textContent = 'Puntos: ' + n + (n && (n < 8 || n > 15) ? ' (recomendado 8–15)' : '');
  $('sPoints').textContent = n; $('mPoints').textContent = n;

  // Estado
  let st = '<span class="badge off">Sin calcular</span>';
  if (m) st = '<span class="badge on">Calculado ✓</span>';
  else if (state.splineOn) st = '<span class="badge bad">No válido</span>';
  $('sState').innerHTML = st;
  $('sSegs').textContent = m ? m.n + ' / ' + m.n : '—';
  $('mIntervals').textContent = m ? m.n : (n >= 2 ? (n - 1) + ' (sin calcular)' : '—');
  $('sLen').textContent = m ? fmt(m.t[m.n], 2) + ' px' : (n >= 2 ? fmt(SplineMath.chordParameters(P)[n - 1], 2) + ' px' : '—');

  // Tabla de puntos
  if (!n) {
    $('ptsBody').innerHTML = '<tr><td colspan="4" style="text-align:center;color:#94a3b8">Sin puntos</td></tr>';
  } else {
    const t = SplineMath.chordParameters(P);
    $('ptsBody').innerHTML = P.map((p, i) =>
      '<tr data-i="' + i + '"' + (i === state.selected ? ' class="sel"' : '') + '><td>' + (i + 1) + '</td><td>' + p.x.toFixed(2) + '</td><td>' + p.y.toFixed(2) + '</td><td>' + t[i].toFixed(2) + '</td></tr>').join('');
  }

  // Tabla de coeficientes
  if (!m) {
    $('coefBody').innerHTML = '<tr><td colspan="11" style="text-align:center;color:#94a3b8">Presiona “CALCULAR SPLINE” para obtener los coeficientes.</td></tr>';
  } else {
    let h = '';
    for (let i = 0; i < m.n; i++) {
      h += '<tr><td>' + i + '</td><td>[' + fmt(m.t[i], 2) + ' , ' + fmt(m.t[i + 1], 2) + ']</td><td>' + fmt(m.h[i], 3) + '</td>' +
        '<td>' + fmt(m.x.a[i]) + '</td><td>' + fmt(m.x.b[i]) + '</td><td>' + fmt(m.x.c[i]) + '</td><td>' + fmt(m.x.d[i]) + '</td>' +
        '<td>' + fmt(m.y.a[i]) + '</td><td>' + fmt(m.y.b[i]) + '</td><td>' + fmt(m.y.c[i]) + '</td><td>' + fmt(m.y.d[i]) + '</td></tr>';
    }
    $('coefBody').innerHTML = h;
  }

  // Verificación
  const v = state.verify;
  if (!v) {
    $('verifBody').innerHTML = '<tr><td style="color:#94a3b8">Calcula el spline para ver los errores.</td></tr>';
  } else {
    const row = (nm, val) => '<tr><td>' + nm + '</td><td>' + val.toExponential(2) + ' px</td><td>' + (val < 1e-6 ? '✓' : '✗') + '</td></tr>';
    $('verifBody').innerHTML =
      row('Error de interpolación máx.', v.interpolacion) + row('Salto C⁰ máx. (posición)', v.C0) +
      row('Salto C¹ máx. (tangente)', v.C1) + row('Salto C² máx. (curvatura)', v.C2) +
      row('Condición natural (X″=Y″=0)', v.frontera_natural);
  }
}
$('ptsBody').addEventListener('mouseover', e => {
  const tr = e.target.closest('tr[data-i]'); if (!tr) return;
  const i = parseInt(tr.dataset.i, 10);
  if (i !== state.hover && state.mode === 'edit') { state.hover = i; scheduleRender(); }
});
$('ptsBody').addEventListener('mouseleave', () => { if (state.hover !== -1) { state.hover = -1; scheduleRender(); } });
$('ptsBody').addEventListener('click', e => {
  const tr = e.target.closest('tr[data-i]'); if (!tr || state.mode !== 'edit') return;
  state.selected = parseInt(tr.dataset.i, 10); scheduleRender(true);
});

/* =====================================================================
   10. INICIO
   ===================================================================== */
if (window.ResizeObserver) new ResizeObserver(resizeCanvas).observe(stage);
window.addEventListener('resize', resizeCanvas);
resizeCanvas();
setMode('edit');
updatePanels();
msg('info', 'Carga una fotografía de la pieza para comenzar.');
