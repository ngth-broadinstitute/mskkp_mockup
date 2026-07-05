const DATA_ROOT = "data";
const PAIR_COLORS = { left: "#2f6f73", right: "#c06938" };
const state = {
  index: null,
  datasets: new Map(),
  geneIndex: [],
  currentGenePayload: null,
  currentCellTypePayload: null,
  cache: new Map(),
  left: null,
  right: null,
};

const els = {};

document.addEventListener("DOMContentLoaded", init);

async function init() {
  cacheElements();
  try {
    const [index, geneIndex] = await Promise.all([
      fetchJson(`${DATA_ROOT}/index.json`),
      fetchJson(`${DATA_ROOT}/genes/index.json`),
    ]);
    state.index = index;
    state.datasets = new Map(index.datasets.map((dataset) => [dataset.id, dataset]));
    state.geneIndex = geneIndex.items || [];
    state.left = index.default_left_dataset;
    state.right = index.default_right_dataset;
    setupControls();
    renderLegend();
    await updateAll();
  } catch (error) {
    renderFatal(error);
  }
}

function cacheElements() {
  Object.assign(els, {
    leftDataset: document.getElementById("left-dataset"),
    rightDataset: document.getElementById("right-dataset"),
    leftTitle: document.getElementById("left-title"),
    rightTitle: document.getElementById("right-title"),
    leftCount: document.getElementById("left-count"),
    rightCount: document.getElementById("right-count"),
    leftUmap: document.getElementById("left-umap"),
    rightUmap: document.getElementById("right-umap"),
    legend: document.getElementById("legend"),
    geneSelector: document.getElementById("gene-selector"),
    geneOptions: document.getElementById("gene-options"),
    geneStatus: document.getElementById("gene-status"),
    genePlot: document.getElementById("gene-plot"),
    leftGeneTable: document.getElementById("left-gene-table"),
    rightGeneTable: document.getElementById("right-gene-table"),
    leftGeneTitle: document.getElementById("left-gene-table-title"),
    rightGeneTitle: document.getElementById("right-gene-table-title"),
    cellTypeSelector: document.getElementById("celltype-selector"),
    cellTypeStatus: document.getElementById("celltype-status"),
    cellTypePlot: document.getElementById("celltype-plot"),
    leftCellTypeTable: document.getElementById("left-celltype-table"),
    rightCellTypeTable: document.getElementById("right-celltype-table"),
    leftCellTypeTitle: document.getElementById("left-celltype-table-title"),
    rightCellTypeTitle: document.getElementById("right-celltype-table-title"),
    leftKey: document.getElementById("left-key"),
    rightKey: document.getElementById("right-key"),
    leftKeyCellType: document.getElementById("left-key-celltype"),
    rightKeyCellType: document.getElementById("right-key-celltype"),
  });
}

function setupControls() {
  const options = state.index.datasets
    .map((dataset) => `<option value="${escapeHtml(dataset.id)}">${escapeHtml(shortDatasetLabel(dataset))}</option>`)
    .join("");
  els.leftDataset.innerHTML = options;
  els.rightDataset.innerHTML = options;
  els.leftDataset.value = state.left;
  els.rightDataset.value = state.right;
  els.leftDataset.addEventListener("change", () => {
    state.left = els.leftDataset.value;
    updateAll();
  });
  els.rightDataset.addEventListener("change", () => {
    state.right = els.rightDataset.value;
    updateAll();
  });

  els.geneOptions.innerHTML = state.geneIndex
    .map((item) => `<option value="${escapeHtml(item.gene)}"></option>`)
    .join("");
  els.geneSelector.value = state.geneIndex[0]?.gene || "CXCL12";
  els.geneSelector.addEventListener("change", updateGeneView);
  els.geneSelector.addEventListener("keydown", (event) => {
    if (event.key === "Enter") updateGeneView();
  });

  els.cellTypeSelector.innerHTML = state.index.cell_types
    .map((cellType) => `<option value="${escapeHtml(cellType)}">${formatLabel(cellType)}</option>`)
    .join("");
  els.cellTypeSelector.value = state.index.cell_types.includes("osteoblast") ? "osteoblast" : state.index.cell_types[0];
  els.cellTypeSelector.addEventListener("change", updateCellTypeView);
}

async function updateAll() {
  const left = state.datasets.get(state.left);
  const right = state.datasets.get(state.right);
  updatePairLabels(left, right);
  await Promise.all([updateOverview(), updateGeneView(), updateCellTypeView()]);
}

function updatePairLabels(left, right) {
  const leftLabel = shortDatasetLabel(left);
  const rightLabel = shortDatasetLabel(right);
  els.leftTitle.textContent = leftLabel;
  els.rightTitle.textContent = rightLabel;
  els.leftKey.textContent = leftLabel;
  els.rightKey.textContent = rightLabel;
  els.leftKeyCellType.textContent = leftLabel;
  els.rightKeyCellType.textContent = rightLabel;
  els.leftGeneTitle.textContent = `${leftLabel} summaries`;
  els.rightGeneTitle.textContent = `${rightLabel} summaries`;
  els.leftCellTypeTitle.textContent = `${leftLabel} summaries`;
  els.rightCellTypeTitle.textContent = `${rightLabel} summaries`;
}

async function updateOverview() {
  const [leftAtlas, rightAtlas] = await Promise.all([
    cachedJson(`${DATA_ROOT}/datasets/${state.left}/atlas.json`),
    cachedJson(`${DATA_ROOT}/datasets/${state.right}/atlas.json`),
  ]);
  renderUmap(els.leftUmap, leftAtlas.cells || []);
  renderUmap(els.rightUmap, rightAtlas.cells || []);
  els.leftCount.textContent = `${formatInteger(leftAtlas.n_cells_total)} cells`;
  els.rightCount.textContent = `${formatInteger(rightAtlas.n_cells_total)} cells`;
}

async function updateGeneView() {
  const gene = resolveGene(els.geneSelector.value.trim());
  if (!gene) return setStatus(els.geneStatus, "Choose a listed gene.");
  els.geneSelector.value = gene.gene;
  setStatus(els.geneStatus, "Loading");
  try {
    const payload = await cachedJson(gene.file);
    state.currentGenePayload = payload;
    const categories = pairedGeneCategories(payload, state.left, state.right);
    renderViolinPlot(els.genePlot, els.geneStatus, categories, "Expression");
    renderSummaryTable(els.leftGeneTable, categories.map((row) => ({ label: row.label, summary: row.leftSummary })), "Cell type");
    renderSummaryTable(els.rightGeneTable, categories.map((row) => ({ label: row.label, summary: row.rightSummary })), "Cell type");
  } catch (error) {
    clearCanvas(els.genePlot);
    setStatus(els.geneStatus, error.message);
  }
}

async function updateCellTypeView() {
  const cellType = els.cellTypeSelector.value;
  setStatus(els.cellTypeStatus, "Loading");
  try {
    const payload = await cachedJson(`${DATA_ROOT}/celltypes/${cleanFilename(cellType)}.json`);
    const genes = (payload.genes || []).slice(0, 12);
    const genePayloads = await Promise.all(genes.map((gene) => cachedJson(`${DATA_ROOT}/genes/${cleanFilename(gene)}.json`)));
    const categories = genePayloads.map((genePayload) => pairedCellTypeCategory(genePayload, cellType, state.left, state.right)).filter(Boolean);
    state.currentCellTypePayload = { cellType, categories };
    renderViolinPlot(els.cellTypePlot, els.cellTypeStatus, categories, "Expression");
    renderSummaryTable(els.leftCellTypeTable, categories.map((row) => ({ label: row.label, summary: row.leftSummary })), "Gene");
    renderSummaryTable(els.rightCellTypeTable, categories.map((row) => ({ label: row.label, summary: row.rightSummary })), "Gene");
  } catch (error) {
    clearCanvas(els.cellTypePlot);
    setStatus(els.cellTypeStatus, error.message);
  }
}

function pairedGeneCategories(payload, leftId, rightId) {
  const leftRows = new Map((payload.datasets[leftId] || []).map((row) => [row.cell_type, row]));
  const rightRows = new Map((payload.datasets[rightId] || []).map((row) => [row.cell_type, row]));
  const order = state.index.cell_types.filter((cellType) => leftRows.has(cellType) || rightRows.has(cellType));
  return order.map((cellType) => ({
    label: formatPlainLabel(cellType),
    leftValues: leftRows.get(cellType)?.values || [],
    rightValues: rightRows.get(cellType)?.values || [],
    leftSummary: leftRows.get(cellType)?.summary || emptySummary(),
    rightSummary: rightRows.get(cellType)?.summary || emptySummary(),
  }));
}

function pairedCellTypeCategory(genePayload, cellType, leftId, rightId) {
  const leftRow = (genePayload.datasets[leftId] || []).find((row) => row.cell_type === cellType);
  const rightRow = (genePayload.datasets[rightId] || []).find((row) => row.cell_type === cellType);
  if (!leftRow && !rightRow) return null;
  return {
    label: genePayload.gene,
    leftValues: leftRow?.values || [],
    rightValues: rightRow?.values || [],
    leftSummary: leftRow?.summary || emptySummary(),
    rightSummary: rightRow?.summary || emptySummary(),
  };
}

function renderSummaryTable(container, rows, firstLabel) {
  if (!rows.length) {
    container.innerHTML = `<div class="empty-state">No rows available.</div>`;
    return;
  }
  container.innerHTML = `
    <table>
      <thead><tr><th>${firstLabel}</th><th class="numeric">Avg</th><th class="numeric">% expr</th><th class="numeric">n</th></tr></thead>
      <tbody>
        ${rows.map((row) => `
          <tr>
            <td>${escapeHtml(row.label)}</td>
            <td class="numeric">${formatNumber(row.summary.avg_expression)}</td>
            <td class="numeric">${formatPercent(row.summary.pct_expressing)}</td>
            <td class="numeric">${formatInteger(row.summary.n)}</td>
          </tr>
        `).join("")}
      </tbody>
    </table>
  `;
}

function renderLegend() {
  els.legend.innerHTML = state.index.cell_types.map((cellType) => `
    <span class="legend-item">
      <span class="swatch" style="background:${state.index.cell_type_colors[cellType] || "#999"}"></span>
      ${formatLabel(cellType)}
    </span>
  `).join("");
}

function renderUmap(canvas, cells) {
  const ctx = prepareCanvas(canvas);
  const width = canvas.clientWidth;
  const height = canvas.clientHeight;
  ctx.clearRect(0, 0, width, height);
  ctx.fillStyle = "#fbfcfc";
  ctx.fillRect(0, 0, width, height);
  if (!cells.length) return;
  const xs = cells.map((cell) => cell.x);
  const ys = cells.map((cell) => cell.y);
  const bounds = { minX: Math.min(...xs), maxX: Math.max(...xs), minY: Math.min(...ys), maxY: Math.max(...ys) };
  const pad = 28;
  cells.forEach((cell) => {
    const x = scale(cell.x, bounds.minX, bounds.maxX, pad, width - pad);
    const y = scale(cell.y, bounds.minY, bounds.maxY, height - pad, pad);
    ctx.globalAlpha = 0.86;
    ctx.fillStyle = state.index.cell_type_colors[cell.cell_type] || "#9aa5a0";
    ctx.beginPath();
    ctx.arc(x, y, 2.45, 0, Math.PI * 2);
    ctx.fill();
  });
  ctx.globalAlpha = 1;
}

function renderViolinPlot(canvas, statusEl, categories, yLabel) {
  const available = categories.filter((row) => row.leftValues.length || row.rightValues.length);
  if (!available.length) {
    clearCanvas(canvas);
    return setStatus(statusEl, "No plot values available.");
  }
  const ctx = prepareCanvas(canvas);
  const width = canvas.clientWidth;
  const height = canvas.clientHeight;
  const margins = { top: 20, right: 18, bottom: 92, left: 58 };
  const plotWidth = width - margins.left - margins.right;
  const plotHeight = height - margins.top - margins.bottom;
  const allValues = available.flatMap((row) => [...row.leftValues, ...row.rightValues]).map(Number).filter(Number.isFinite);
  const yMax = niceCeil(Math.max(0.05, quantile(allValues, 0.995)));
  ctx.clearRect(0, 0, width, height);
  ctx.fillStyle = "#fbfcfc";
  ctx.fillRect(0, 0, width, height);
  drawAxes(ctx, margins, plotWidth, plotHeight, yMax, yLabel);
  const band = plotWidth / available.length;
  const groupWidth = Math.min(58, band * 0.68);
  const offset = Math.min(20, groupWidth * 0.28);
  available.forEach((row, index) => {
    const center = margins.left + band * (index + 0.5);
    drawDistribution(ctx, row.leftValues, center - offset, groupWidth * 0.38, yMax, margins, plotHeight, PAIR_COLORS.left);
    drawDistribution(ctx, row.rightValues, center + offset, groupWidth * 0.38, yMax, margins, plotHeight, PAIR_COLORS.right);
    drawCategoryLabel(ctx, row.label, center, margins.top + plotHeight + 14);
  });
  setStatus(statusEl, "");
}

function drawAxes(ctx, margins, plotWidth, plotHeight, yMax, yLabel) {
  ctx.strokeStyle = "#c4d0ca";
  ctx.lineWidth = 1;
  ctx.beginPath();
  ctx.moveTo(margins.left, margins.top);
  ctx.lineTo(margins.left, margins.top + plotHeight);
  ctx.lineTo(margins.left + plotWidth, margins.top + plotHeight);
  ctx.stroke();
  ctx.fillStyle = "#63706b";
  ctx.font = "12px Inter, sans-serif";
  ctx.textAlign = "right";
  ctx.textBaseline = "middle";
  for (let i = 0; i <= 4; i += 1) {
    const value = (yMax * i) / 4;
    const y = yScale(value, yMax, margins.top + plotHeight, margins.top);
    ctx.strokeStyle = i === 0 ? "#c4d0ca" : "#e5ece8";
    ctx.beginPath();
    ctx.moveTo(margins.left, y);
    ctx.lineTo(margins.left + plotWidth, y);
    ctx.stroke();
    ctx.fillText(formatAxisNumber(value), margins.left - 8, y);
  }
  ctx.save();
  ctx.translate(16, margins.top + plotHeight / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.textAlign = "center";
  ctx.fillText(yLabel, 0, 0);
  ctx.restore();
}

function drawDistribution(ctx, rawValues, x, maxHalfWidth, yMax, margins, plotHeight, color) {
  const values = rawValues.map(Number).filter(Number.isFinite).sort((a, b) => a - b);
  if (!values.length) return;
  if (values.length < 8 || new Set(values.map((v) => v.toFixed(4))).size < 3) {
    drawJitter(ctx, values, x, maxHalfWidth, yMax, margins, plotHeight, color);
    return drawMedian(ctx, values, x, yMax, margins, plotHeight, color);
  }
  const bins = 34;
  const bandwidth = Math.max(yMax / 18, 0.04);
  const points = [];
  let maxDensity = 0;
  for (let i = 0; i < bins; i += 1) {
    const value = (yMax * i) / (bins - 1);
    const density = values.reduce((sum, sample) => {
      const z = (value - sample) / bandwidth;
      return sum + Math.exp(-0.5 * z * z);
    }, 0) / values.length;
    maxDensity = Math.max(maxDensity, density);
    points.push({ value, density });
  }
  ctx.beginPath();
  points.forEach((point, index) => {
    const w = (point.density / maxDensity) * maxHalfWidth;
    const y = yScale(point.value, yMax, margins.top + plotHeight, margins.top);
    if (index === 0) ctx.moveTo(x - w, y);
    else ctx.lineTo(x - w, y);
  });
  [...points].reverse().forEach((point) => {
    const w = (point.density / maxDensity) * maxHalfWidth;
    const y = yScale(point.value, yMax, margins.top + plotHeight, margins.top);
    ctx.lineTo(x + w, y);
  });
  ctx.closePath();
  ctx.fillStyle = hexToRgba(color, 0.32);
  ctx.strokeStyle = color;
  ctx.lineWidth = 1.2;
  ctx.fill();
  ctx.stroke();
  drawMedian(ctx, values, x, yMax, margins, plotHeight, color);
}

function drawJitter(ctx, values, x, maxHalfWidth, yMax, margins, plotHeight, color) {
  ctx.fillStyle = hexToRgba(color, 0.62);
  values.slice(0, 80).forEach((value, index) => {
    const jitter = (Math.sin((index + 1) * 12.9898) * 43758.5453 % 1) * maxHalfWidth * 1.5;
    const y = yScale(value, yMax, margins.top + plotHeight, margins.top);
    ctx.beginPath();
    ctx.arc(x + jitter, y, 2.2, 0, Math.PI * 2);
    ctx.fill();
  });
}

function drawMedian(ctx, values, x, yMax, margins, plotHeight, color) {
  const y = yScale(quantile(values, 0.5), yMax, margins.top + plotHeight, margins.top);
  ctx.fillStyle = "#fff";
  ctx.strokeStyle = color;
  ctx.lineWidth = 2;
  ctx.beginPath();
  ctx.arc(x, y, 3.5, 0, Math.PI * 2);
  ctx.fill();
  ctx.stroke();
}

function drawCategoryLabel(ctx, label, x, y) {
  ctx.save();
  ctx.fillStyle = "#43504b";
  ctx.font = "11px Inter, sans-serif";
  ctx.textAlign = "right";
  ctx.textBaseline = "top";
  ctx.translate(x - 4, y + 50);
  ctx.rotate(-Math.PI / 4);
  ctx.fillText(truncate(label, 24), 0, 0);
  ctx.restore();
}

async function fetchJson(url) {
  const response = await fetch(url);
  if (!response.ok) throw new Error(`Could not load ${url}`);
  return response.json();
}

async function cachedJson(url) {
  if (!state.cache.has(url)) state.cache.set(url, fetchJson(url));
  return state.cache.get(url);
}

function resolveGene(value) {
  return state.geneIndex.find((item) => item.gene.toLowerCase() === value.toLowerCase()) || state.geneIndex[0];
}

function shortDatasetLabel(dataset) {
  return `${dataset.geoSeries} ${dataset.species === "Homo sapiens" ? "Human" : "Mouse"}`;
}

function prepareCanvas(canvas) {
  const dpr = window.devicePixelRatio || 1;
  const rect = canvas.getBoundingClientRect();
  const width = Math.max(1, Math.round(rect.width));
  const height = Math.max(1, Math.round(rect.height));
  if (canvas.width !== width * dpr || canvas.height !== height * dpr) {
    canvas.width = width * dpr;
    canvas.height = height * dpr;
  }
  const ctx = canvas.getContext("2d");
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  return ctx;
}

function clearCanvas(canvas) {
  const ctx = prepareCanvas(canvas);
  ctx.clearRect(0, 0, canvas.clientWidth, canvas.clientHeight);
}

function renderFatal(error) {
  document.body.innerHTML = `<main class="page-shell"><div class="section-block"><div class="empty-state">${escapeHtml(error.message)}</div></div></main>`;
}

function setStatus(el, message) { el.textContent = message || ""; }
function emptySummary() { return { n: 0, avg_expression: 0, pct_expressing: 0, median: 0 }; }
function scale(value, d0, d1, r0, r1) { return d1 === d0 ? (r0 + r1) / 2 : r0 + ((value - d0) / (d1 - d0)) * (r1 - r0); }
function yScale(value, yMax, r0, r1) { return scale(Math.max(0, Math.min(yMax, value || 0)), 0, yMax, r0, r1); }
function quantile(values, q) {
  const sorted = values.map(Number).filter(Number.isFinite).sort((a, b) => a - b);
  if (!sorted.length) return 0;
  const pos = (sorted.length - 1) * q;
  const base = Math.floor(pos);
  const rest = pos - base;
  return sorted[base + 1] === undefined ? sorted[base] : sorted[base] + rest * (sorted[base + 1] - sorted[base]);
}
function niceCeil(value) {
  const exponent = Math.floor(Math.log10(value));
  const fraction = value / 10 ** exponent;
  const nice = fraction <= 1 ? 1 : fraction <= 2 ? 2 : fraction <= 5 ? 5 : 10;
  return nice * 10 ** exponent;
}
function hexToRgba(hex, alpha) {
  const clean = hex.replace("#", "");
  return `rgba(${parseInt(clean.slice(0, 2), 16)}, ${parseInt(clean.slice(2, 4), 16)}, ${parseInt(clean.slice(4, 6), 16)}, ${alpha})`;
}
function formatLabel(value) { return escapeHtml(formatPlainLabel(value)); }
function formatPlainLabel(value) { return String(value || "").replace(/_/g, " "); }
function formatInteger(value) { return new Intl.NumberFormat("en-US").format(Number(value) || 0); }
function formatNumber(value) { return Number(value || 0).toFixed(3); }
function formatPercent(value) { return `${Math.round(Number(value || 0) * 100)}%`; }
function formatAxisNumber(value) { return value >= 10 ? String(Math.round(value)) : value >= 1 ? value.toFixed(1) : value.toFixed(2); }
function truncate(value, max) { return String(value || "").length > max ? `${String(value).slice(0, max - 1)}...` : String(value || ""); }
function cleanFilename(value) { return String(value).replace(/[^A-Za-z0-9_.-]+/g, "_"); }
function escapeHtml(value) {
  return String(value).replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;").replace(/"/g, "&quot;").replace(/'/g, "&#039;");
}

let resizeFrame = null;
window.addEventListener("resize", () => {
  window.cancelAnimationFrame(resizeFrame);
  resizeFrame = window.requestAnimationFrame(() => {
    updateOverview();
    if (state.currentGenePayload) {
      renderViolinPlot(els.genePlot, els.geneStatus, pairedGeneCategories(state.currentGenePayload, state.left, state.right), "Expression");
    }
    if (state.currentCellTypePayload) {
      renderViolinPlot(els.cellTypePlot, els.cellTypeStatus, state.currentCellTypePayload.categories, "Expression");
    }
  });
});
