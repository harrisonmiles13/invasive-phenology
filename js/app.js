/* Invasive Phenology Atlas: maps and a location panel built from site/data/
 * (see pipeline/build_site.py for the file formats). */
(function () {
  "use strict";
  const M = window.PhenoModel;
  const DATA = "data/";
  const EXT = window.PHENO_DATA_EXT || ".bin"; // hosts without binary files get base64 text
  const $ = (s, el = document) => el.querySelector(s);
  const h = (tag, attrs = {}, ...kids) => {
    const el = document.createElement(tag);
    for (const [k, v] of Object.entries(attrs)) {
      if (k === "class") el.className = v;
      else if (k.startsWith("on")) el.addEventListener(k.slice(2), v);
      else if (v !== false && v != null) el.setAttribute(k, v === true ? "" : v);
    }
    for (const kid of kids.flat()) if (kid != null && kid !== false) el.append(kid.nodeType ? kid : document.createTextNode(kid));
    return el;
  };

  // ---------------------------------------------------------------- state
  const state = {
    mode: "today", driver: "agdd50", tview: "value",
    sp: "slf", ev: null, year: "current", view: "date", units: "F",
    rangeOnly: true, showObs: true, stageDoy: null, stageSrc: "current",
    loc: { lat: 40.42, lon: -75.93, name: "Reading, Pennsylvania (example)" },
  };
  let man = null, validation = {};

  // ---------------------------------------------------------------- data
  const cache = new Map();
  function fetchBin(path) {
    const url = DATA + path.replace(/\.bin$/, EXT);
    if (!cache.has(url)) {
      cache.set(url, fetch(url).then(async (r) => {
        if (!r.ok) throw new Error(`${path}: HTTP ${r.status}`);
        let buf = EXT.endsWith(".txt")
          ? Uint8Array.from(atob((await r.text()).trim()), (c) => c.charCodeAt(0)).buffer
          : await r.arrayBuffer();
        const b = new Uint8Array(buf, 0, 2);
        if (b[0] === 0x1f && b[1] === 0x8b) buf = await new Response(new Blob([buf]).stream().pipeThrough(new DecompressionStream("gzip"))).arrayBuffer();
        return new Uint16Array(buf);
      }).catch((e) => { cache.delete(url); throw e; }));
    }
    return cache.get(url);
  }
  const fetchJSON = (p) => fetch(DATA + p).then((r) => { if (!r.ok) throw new Error(`${p}: HTTP ${r.status}`); return r.json(); });
  const N = () => man.grid.nrows * man.grid.ncols;
  const spById = (id) => man.species.find((s) => s.id === id);
  const evOf = (sp, id) => sp.events.find((e) => e.id === id);
  const cur = () => man.current;
  const curYear = () => cur().year;
  const todayDoy = () => cur().today_doy;
  const ev = (part) => `ev/${state.sp}/${state.ev}.${part}.bin`;

  function cellOf(lat, lon) {
    const g = man.grid;
    const i = Math.floor((g.north - lat) / g.dlat), j = Math.floor((lon - g.west) / g.dlon);
    if (i < 0 || j < 0 || i >= g.nrows || j >= g.ncols) return null;
    return { i, j, idx: i * g.ncols + j };
  }

  // ---------------------------------------------------------------- units
  const ddUnit = () => (state.units === "F" ? "°F·days" : "°C·days");
  const toUnit = (c) => (state.units === "F" ? c * 1.8 : c);
  const fmtDD = (c) => Math.round(toUnit(c)).toLocaleString();
  const yearLabel = (y) => (y === "current" ? String(curYear()) : String(y));

  // ---------------------------------------------------------------- colour
  const hex = (s) => [1, 3, 5].map((k) => parseInt(s.slice(k, k + 2), 16));
  const VIRIDIS = ["#440154", "#482878", "#3e4989", "#31688e", "#26828e", "#1f9e89", "#35b779", "#6ece58", "#b5de2b", "#fde725"].map(hex);
  const ramp = (stops, t) => {
    const x = Math.min(1, Math.max(0, t)) * (stops.length - 1), k = Math.min(stops.length - 2, Math.floor(x)), f = x - k;
    return stops[k].map((c, n) => Math.round(c + (stops[k + 1][n] - c) * f));
  };
  const DIVERGE = ["#a50f15", "#de2d26", "#fc9272", "#fdd5c4", "#e6e4df", "#d1e5f0", "#92c5de", "#4393c3", "#1f5f9f"].map(hex);
  const ORANGES = ["#fee6c6", "#fdc58c", "#fd9f55", "#f2701d", "#d4500b", "#a33a05", "#6e2a05"].map(hex);
  const BLUES = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"].map(hex);
  const HEAT = ["#fff3c4", "#fee391", "#fec44f", "#fe9929", "#ec7014", "#cc4c02", "#993404", "#662506"].map(hex);
  const GREY = hex("#bdb9b0"), PALE = hex("#e4e2dc");
  const SPECIAL = {
    notyet: { c: PALE, label: "Not expected this year" },
    notreached: { c: GREY, label: "Not reached" },
    nonormal: { c: GREY, label: "Reached in < half of years" },
    beyond: { c: hex("#67000d"), label: "Past a normal full year's heat" },
  };
  const PATTERN = { 3: "dots", 4: "hatch" }; // forecast / projected cells

  // Calendar bins: half-months or months spanning the 2nd-98th percentile.
  function dateBins(values) {
    const v = values.filter((x) => x > 0).sort((a, b) => a - b);
    if (!v.length) return [91, 121];
    const p = (q) => v[Math.floor(q * (v.length - 1))];
    const lo = p(0.02), hi = p(0.98), semi = [], month = [];
    for (let m = 0; m < 12; m++) {
      semi.push(M.dateToDoy(new Date(Date.UTC(2021, m, 1))), M.dateToDoy(new Date(Date.UTC(2021, m, 16))));
      month.push(M.dateToDoy(new Date(Date.UTC(2021, m, 1))));
    }
    semi.push(366); month.push(366);
    for (const edges of [semi, month]) {
      const a = edges.filter((e) => e <= lo).pop() ?? edges[0], b = edges.find((e) => e > hi) ?? 366;
      const sel = edges.filter((e) => e >= a && e <= b);
      if (sel.length - 1 <= 12 || edges === month) return sel;
    }
  }
  function dateScale(edges) {
    const n = edges.length - 1, classes = [];
    for (let k = 0; k < n; k++) {
      const a = M.fmtDoy(edges[k]), b = M.fmtDoy(edges[k + 1] - 1);
      classes.push({ c: ramp(VIRIDIS, 1 - k / Math.max(1, n - 1)), label: k === 0 ? `by ${b}` : k === n - 1 ? `${a} or later` : `${a}–${b}` });
    }
    return { classes, cls: (d) => { for (let k = 1; k < n; k++) if (d < edges[k]) return k - 1; return n - 1; } };
  }
  function stepScale(edges, colors, fmt) {
    const classes = [];
    for (let k = 0; k <= edges.length; k++) {
      const a = edges[k - 1], b = edges[k];
      classes.push({ c: colors[k], label: k === 0 ? `< ${fmt(b)}` : k === edges.length ? `${fmt(a)}+` : `${fmt(a)}–${fmt(b)}` });
    }
    return { classes, cls: (v) => { let k = 0; while (k < edges.length && v >= edges[k]) k++; return k; } };
  }
  const niceStep = (span, n) => { const raw = span / n, p = 10 ** Math.floor(Math.log10(raw)); return [1, 2, 2.5, 5, 10].map((m) => m * p).find((s) => s >= raw); };

  // ---------------------------------------------------------------- views
  const SPECIES_VIEWS = {
    date: "Date", anom: "vs. normal", normal: "Normal", sd: "Variability", pct: "Years reached", status: "Status on date",
  };
  const TODAY_VIEWS = { value: "Accumulated", anom: "vs. normal", ahead: "Days ahead / behind" };
  const STAGE_ACTIONS = {
    "Breaking Leaf Buds": "Leaf-out. Invasive shrubs often leaf out before the native canopy, when they are easiest to find and treat.",
    "Flower Budding": "Buds forming. Plan removal before flowering and seed set.",
    Flowering: "Flowering. Plants are easiest to identify; cut, pull or treat before seed matures.",
    Fruiting: "Fruit and seed forming. Bag seed heads when removing plants so seed is not spread.",
    "Colored Leaves": "Fall color. Invasive shrubs often hold leaves after natives drop them: a fall foliar treatment window.",
    Larva: "Larvae feeding. Damage appears and treatments that target young larvae work best early.",
    Nymph: "Nymphs active. Monitor with traps; young nymphs are the easiest stage to control.",
    Adult: "Adults active. Trapping, dispersal and egg laying.",
    Egg: "Eggs present. Survey for and destroy egg masses.",
  };

  let landMask = null;
  const isLand = (i) => landMask && landMask[i] === 1;
  let rangeMask = null;

  async function buildView() {
    if (state.mode === "today") return buildToday();
    const sp = spById(state.sp), e = evOf(sp, state.ev);
    const n = N(), np = man.normal_period.join("–");
    const isCur = state.year === "current";
    const nrm = await fetchBin(ev("normal"));
    const nAt = (i) => nrm[i] / 10;
    const label = `${sp.common}: ${e.label.toLowerCase()}`;
    if (state.view === "status") return buildStatus(sp);
    if (state.view === "normal" || state.view === "date") {
      const sc = dateScale(dateBins(Array.from(nrm.subarray(0, n), (v) => v / 10)));
      if (state.view === "normal") return {
        title: `${label}, normal date`, sub: `Mean date, ${np}`, scale: sc, specials: ["nonormal"],
        at(i) {
          if (!isLand(i)) return null;
          const m = nAt(i);
          if (!m) return { k: "nonormal", text: `reached in ${nrm[2 * n + i]}% of ${np} years` };
          return { k: sc.cls(m), text: `normally ${M.fmtDoy(m)} (± ${(nrm[n + i] / 10).toFixed(0)} d)` };
        },
      };
      const a = await fetchBin(ev(isCur ? "current" : String(state.year)));
      return {
        title: `${label}, ${yearLabel(state.year)}`,
        sub: isCur ? `Observed through ${fmtISO(cur().obs_through)}; NWS forecast to ${M.fmtDoy(cur().forecast_through_doy, curYear())}; normal temperatures after that` : "Date reached",
        scale: sc, specials: [isCur ? "notyet" : "notreached"], patterns: isCur,
        at(i) {
          if (!isLand(i)) return null;
          const { doy, src } = isCur ? M.decodeCurrent(a[i]) : { doy: a[i], src: 1 };
          if (!doy) return { k: isCur ? "notyet" : "notreached", text: isCur ? "not expected this year" : "not reached" };
          return { k: sc.cls(doy), pat: PATTERN[src], text: `${M.fmtDoy(doy, yearNum())} (${srcWord(src)})` };
        },
      };
    }
    if (state.view === "anom") {
      const a = await fetchBin(ev(isCur ? "current" : String(state.year)));
      const edges = [-21, -14, -7, -3, 3, 7, 14, 21];
      const labels = ["21+ days early", "14–21 early", "7–14 early", "3–7 early", "Within 3 days", "3–7 late", "7–14 late", "14–21 late", "21+ days late"];
      const classes = labels.map((l, k) => ({ c: DIVERGE[k], label: l }));
      const cls = (d) => { let k = 0; while (k < edges.length && d > edges[k]) k++; return k; };
      return {
        title: `${label}, ${yearLabel(state.year)} vs. normal`, sub: `Days earlier or later than the ${np} mean`,
        scale: { classes, cls }, specials: ["nonormal"], diverging: true, patterns: isCur,
        at(i) {
          if (!isLand(i)) return null;
          const m = nAt(i);
          const { doy, src } = isCur ? M.decodeCurrent(a[i]) : { doy: a[i], src: 1 };
          if (!m) return { k: "nonormal", text: "no normal date here" };
          if (!doy) return { k: classes.length - 1, text: `not expected this year (normally ${M.fmtDoy(m)})` };
          const d = doy - m;
          return { k: cls(d), pat: PATTERN[src], text: `${M.fmtDoy(doy, yearNum())} (${srcWord(src)}), ${Math.abs(d).toFixed(0)} days ${d < 0 ? "early" : "late"}` };
        },
      };
    }
    if (state.view === "sd") {
      const sc = stepScale([3, 5, 7, 10, 14, 21], ORANGES, (v) => `${v} d`);
      return {
        title: `${label}, year-to-year variability`, sub: `Standard deviation of the date, ${np}`, scale: sc, specials: ["nonormal"],
        at(i) { if (!isLand(i)) return null; if (!nrm[i]) return { k: "nonormal", text: "rarely reached" }; const s = nrm[n + i] / 10; return { k: sc.cls(s), text: `± ${s.toFixed(1)} days` }; },
      };
    }
    if (state.view === "pct") {
      const labels = ["Never", "1–24%", "25–49%", "50–74%", "75–89%", "90–99%", "Every year"];
      const classes = labels.map((l, k) => ({ c: k ? BLUES[k] : GREY, label: l }));
      const cls = (p) => (p === 0 ? 0 : p < 25 ? 1 : p < 50 ? 2 : p < 75 ? 3 : p < 90 ? 4 : p < 100 ? 5 : 6);
      return {
        title: `${label}, share of years reached`, sub: `Share of ${np} years in which it happened by Dec 31`, scale: { classes, cls }, specials: [],
        at(i) { if (!isLand(i)) return null; const g = nrm[2 * n + i]; return { k: cls(g), text: `${g}% of years` }; },
      };
    }
  }

  async function buildStatus(sp) {
    const np = man.normal_period.join("–");
    const src = state.stageSrc, d = state.stageDoy;
    const layers = await Promise.all(sp.events.map((e) => fetchBin(`ev/${sp.id}/${e.id}.${src === "normal" ? "normal" : "current"}.bin`)));
    const n = N();
    const val = (k, i) => (src === "normal" ? layers[k][i] / 10 : M.decodeCurrent(layers[k][i]).doy);
    // Order events by their typical date across the map.
    const med = sp.events.map((e, k) => {
      const v = []; for (let i = 0; i < n; i += 37) { const x = val(k, i); if (x > 0) v.push(x); }
      v.sort((a, b) => a - b); return v.length ? v[v.length >> 1] : 999;
    });
    const order = sp.events.map((e, k) => k).sort((a, b) => med[a] - med[b]);
    const classes = [{ c: PALE, label: `Before ${sp.events[order[0]].label.toLowerCase()}` }]
      .concat(order.map((k, j) => ({ c: ramp(VIRIDIS, 1 - j / Math.max(1, order.length - 1)), label: `Since ${sp.events[k].label.toLowerCase()}` })));
    const when = src === "normal" ? `${M.fmtDoy(d)} in a normal year (${np})` : `${M.fmtDoy(d, curYear())}, ${curYear()}`;
    return {
      title: `${sp.common} on ${when}`, sub: "Latest event reached, all modeled events for this species", scale: { classes, cls: (k) => k }, specials: [],
      at(i) {
        if (!isLand(i)) return null;
        let s = 0;
        order.forEach((k, j) => { const v = val(k, i); if (v > 0 && v <= d) s = j + 1; });
        return { k: s, text: classes[s].label.toLowerCase() };
      },
    };
  }

  const generalCache = {};
  async function generalLayer(id) {
    if (!generalCache[id]) generalCache[id] = fetchBin(`general/${id}.bin`);
    return generalCache[id];
  }
  async function buildToday() {
    const g = await generalLayer(state.driver), n = N();
    const drv = man.drivers.find((d) => d.id === state.driver);
    const ok = (i) => g[i] !== 65535 && isLand(i);
    const today = `${M.fmtDoy(todayDoy(), curYear())}, ${curYear()}`;
    if (state.tview === "value") {
      const vals = []; for (let i = 0; i < n; i += 7) if (ok(i)) vals.push(toUnit(g[i]));
      vals.sort((a, b) => a - b);
      const hi = vals[Math.floor(vals.length * 0.98)] || 1, step = niceStep(hi, 9), edges = [];
      for (let x = step; x < hi; x += step) edges.push(x);
      const colors = edges.concat([0]).map((_, k) => ramp(HEAT, k / Math.max(1, edges.length)));
      const sc = stepScale(edges, colors, (v) => Math.round(v).toLocaleString());
      return {
        title: `${drv.label}: accumulated to ${today}`, sub: `Since Jan 1, ${ddUnit()}. Includes today's NWS forecast.`, scale: sc, specials: [],
        at(i) { if (!ok(i)) return null; return { k: sc.cls(toUnit(g[i])), text: `${fmtDD(g[i])} ${ddUnit()} (normal ${fmtDD(g[n + i])})` }; },
      };
    }
    if (state.tview === "anom") {
      const e = [-400, -200, -100, -50, 50, 100, 200, 400].map((x) => (state.units === "F" ? x : Math.round(x / 1.8)));
      const labels = [`${-e[0]}+ below`, `${-e[1]}–${-e[0]} below`, `${-e[2]}–${-e[1]} below`, `${-e[3]}–${-e[2]} below`, "Near normal", `${e[4]}–${e[5]} above`, `${e[5]}–${e[6]} above`, `${e[6]}–${e[7]} above`, `${e[7]}+ above`];
      const classes = labels.map((l, k) => ({ c: DIVERGE[8 - k], label: l }));
      const cls = (d) => { let k = 0; while (k < e.length && d > e[k]) k++; return k; };
      return {
        title: `${drv.label}: ${today} vs. normal`, sub: `${ddUnit()} above or below the ${man.normal_period.join("–")} normal to date`,
        scale: { classes, cls }, specials: [], diverging: true,
        at(i) { if (!ok(i)) return null; const d = toUnit(g[i] - g[n + i]); return { k: cls(d), text: `${d >= 0 ? "+" : "−"}${Math.abs(Math.round(d)).toLocaleString()} ${ddUnit()} vs. normal` }; },
      };
    }
    const edges = [-21, -14, -7, -3, 3, 7, 14, 21];
    const labels = ["21+ days behind", "14–21 behind", "7–14 behind", "3–7 behind", "Within 3 days", "3–7 ahead", "7–14 ahead", "14–21 ahead", "21+ days ahead"];
    const classes = labels.map((l, k) => ({ c: DIVERGE[8 - k], label: l }));
    const cls = (d) => { let k = 0; while (k < edges.length && d > edges[k]) k++; return k; };
    return {
      title: `Season ${today}: days ahead or behind`, sub: `Heat accumulated so far compared with the date the ${man.normal_period.join("–")} normal reaches it (${drv.label.toLowerCase()})`,
      scale: { classes, cls }, specials: ["beyond"], diverging: true,
      at(i) {
        if (!ok(i) || g[2 * n + i] === 65535) return null;
        if (g[2 * n + i] === 3000) return { k: "beyond", text: "already more heat than a normal full year" };
        const d = g[2 * n + i] - 1000;
        return { k: cls(d), text: d === 0 ? "on schedule" : `${Math.abs(d)} days ${d > 0 ? "ahead of" : "behind"} normal` };
      },
    };
  }

  const yearNum = () => (state.year === "current" ? curYear() : +state.year);
  const srcWord = (s) => M.SOURCES[s] || "";
  const fmtISO = (iso) => { const d = new Date(iso + "T00:00:00Z"); return `${M.MONTHS[d.getUTCMonth()]} ${d.getUTCDate()}`; };

  // ---------------------------------------------------------------- raster
  const mercY = (lat) => Math.log(Math.tan(Math.PI / 4 + (lat * Math.PI) / 360));
  const gridBounds = () => { const g = man.grid; return [[g.north - g.nrows * g.dlat, g.west], [g.north, g.west + g.ncols * g.dlon]]; };
  function rasterURL(v) {
    const g = man.grid, [[s], [nth]] = gridBounds(), W = g.ncols * 2, yN = mercY(nth), yS = mercY(s);
    const H = Math.round((W * (yN - yS)) / ((g.ncols * g.dlon * Math.PI) / 180));
    const n = N(), rgba = new Uint8ClampedArray(n * 4), pat = new Uint8Array(n);
    const dim = state.mode === "species" && rangeMask;
    for (let i = 0; i < n; i++) {
      const r = v.at(i);
      if (!r) continue;
      const c = typeof r.k === "string" ? SPECIAL[r.k].c : v.scale.classes[r.k].c;
      let a = 255;
      if (dim && !rangeMask[i]) { if (state.rangeOnly) continue; a = 70; }
      rgba.set([c[0], c[1], c[2], a], i * 4);
      if (r.pat === "dots") pat[i] = 1; else if (r.pat === "hatch") pat[i] = 2;
    }
    const cv = document.createElement("canvas"); cv.width = W; cv.height = H;
    const ctx = cv.getContext("2d"), img = ctx.createImageData(W, H);
    for (let r = 0; r < H; r++) {
      const y = yN - ((r + 0.5) / H) * (yN - yS), lat = (Math.atan(Math.exp(y)) * 360) / Math.PI - 90;
      const base = Math.min(g.nrows - 1, Math.floor((g.north - lat) / g.dlat)) * g.ncols;
      for (let x = 0; x < W; x++) {
        const i = base + (x >> 1), o = i * 4, p = (r * W + x) * 4;
        img.data[p] = rgba[o]; img.data[p + 1] = rgba[o + 1]; img.data[p + 2] = rgba[o + 2]; img.data[p + 3] = rgba[o + 3];
        if (rgba[o + 3] && ((pat[i] === 2 && (x + r) % 5 === 0) || (pat[i] === 1 && x % 3 === 0 && r % 3 === 0))) {
          img.data[p] = img.data[p] * 0.55 + 255 * 0.45; img.data[p + 1] = img.data[p + 1] * 0.55 + 255 * 0.45; img.data[p + 2] = img.data[p + 2] * 0.55 + 255 * 0.45;
        }
      }
    }
    ctx.putImageData(img, 0, 0);
    return cv.toDataURL();
  }

  // ---------------------------------------------------------------- map
  const darkMode = () => { const t = document.documentElement.getAttribute("data-theme"); return t ? t === "dark" : matchMedia("(prefers-color-scheme: dark)").matches; };
  const cssVar = (n) => getComputedStyle(document.documentElement).getPropertyValue(n).trim();
  const tileURL = (k) => `https://{s}.basemaps.cartocdn.com/${darkMode() ? "dark" : "light"}_${k}/{z}/{x}/{y}{r}.png`;
  let map, overlay, base, labels, lines, marker, obsLayer, obsRenderer, currentView = null, boundaries = null;

  function initMap() {
    map = L.map("map", { zoomSnap: 0.25, minZoom: 3, maxZoom: 11, maxBounds: [[5, -175], [72, -35]] });
    for (const [name, z] of [["raster", 350], ["lines", 420], ["labels", 430], ["obs", 440]]) {
      const p = map.createPane(name); p.style.zIndex = z; if (name !== "raster") p.style.pointerEvents = "none";
    }
    map.attributionControl.setPrefix(false);
    base = L.tileLayer(tileURL("nolabels"), { subdomains: "abcd", maxZoom: 19, attribution: '&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> &copy; <a href="https://carto.com/attributions">CARTO</a>' }).addTo(map);
    labels = L.tileLayer(tileURL("only_labels"), { subdomains: "abcd", maxZoom: 19, pane: "labels" }).addTo(map);
    overlay = L.imageOverlay("data:image/gif;base64,R0lGODlhAQABAAAAACw=", gridBounds(), { pane: "raster", opacity: 0.9, interactive: false }).addTo(map);
    if (boundaries) lines = L.polyline(boundaries, { renderer: L.canvas({ pane: "lines" }), color: cssVar("--map-line"), weight: 0.7, interactive: false }).addTo(map);
    obsRenderer = L.canvas({ pane: "obs" });
    obsLayer = L.layerGroup().addTo(map);
    marker = L.circleMarker([state.loc.lat, state.loc.lon], { radius: 7, weight: 3, color: cssVar("--accent"), fillColor: "#fff", fillOpacity: 1, interactive: false, pane: "markerPane" }).addTo(map);
    map.fitBounds([[24.6, -124.5], [49.2, -67]], { padding: [4, 4] });
    map.on("mousemove", (e) => hover(e.latlng));
    map.on("click", (e) => { hover(e.latlng); selectLocation(e.latlng.lat, e.latlng.lng); });
    map.on("zoomend", sizeObs);
  }
  function restyle() {
    base.setUrl(tileURL("nolabels")); labels.setUrl(tileURL("only_labels"));
    if (lines) lines.setStyle({ color: cssVar("--map-line") });
    marker.setStyle({ color: cssVar("--accent") });
  }
  function hover(ll) {
    const c = man && currentView ? cellOf(ll.lat, ll.lng) : null;
    const r = c ? currentView.at(c.idx) : null;
    const where = `${Math.abs(ll.lat).toFixed(2)}°${ll.lat >= 0 ? "N" : "S"} ${Math.abs(ll.lng).toFixed(2)}°W`;
    let extra = "";
    if (r && state.mode === "species" && rangeMask && c && !rangeMask[c.idx]) extra = " · no iNaturalist records nearby";
    $("#readout").textContent = r ? `${where}: ${r.text}${extra}` : `${where}: outside the mapped area (contiguous US)`;
  }

  let renderToken = null;
  async function renderMap() {
    const token = (renderToken = {});
    $("#map-wrap").classList.add("loading");
    try {
      if (state.mode === "species") await loadRange();
      const v = await buildView();
      if (token !== renderToken) return;
      currentView = v;
      $("#map-title").textContent = v.title;
      $("#map-sub").textContent = v.sub;
      overlay.setUrl(rasterURL(v));
      drawLegend(v);
      drawObs();
    } catch (e) {
      $("#map-sub").textContent = `Could not load map data (${e.message}). Reload to try again.`;
    } finally {
      $("#map-wrap").classList.remove("loading");
    }
  }
  function drawLegend(v) {
    const el = $("#legend");
    el.replaceChildren();
    const sw = (c, label, cls = "") => h("li", {}, h("span", { class: "sw " + cls, style: `background-color: rgb(${c.join(",")})` }), h("span", {}, label));
    const list = h("ul", { class: "legend-list" });
    v.scale.classes.forEach((k) => list.append(sw(k.c, k.label)));
    for (const s of v.specials || []) list.append(sw(SPECIAL[s].c, SPECIAL[s].label));
    el.append(list);
    if (v.patterns) el.append(h("ul", { class: "legend-list legend-src" },
      h("li", {}, h("span", { class: "sw solid" }), h("span", {}, "Observed")),
      h("li", {}, h("span", { class: "sw dots" }), h("span", {}, "NWS forecast")),
      h("li", {}, h("span", { class: "sw hatch" }), h("span", {}, "Projected with normal temperatures"))));
    if (state.mode === "species" && rangeMask) el.append(h("p", { class: "note" }, state.rangeOnly
      ? "Shown only within about 25 km of iNaturalist records of this species."
      : "Faded where iNaturalist has no records of this species nearby."));
  }

  async function loadRange() {
    const sp = spById(state.sp);
    rangeMask = null;
    if (!sp.has_range) return;
    try { rangeMask = await fetchBin(`range/${sp.id}.bin`); } catch { rangeMask = null; }
  }

  // iNaturalist records for the selected species and stage.
  const obsCache = {};
  async function drawObs() {
    obsLayer.clearLayers();
    if (state.mode !== "species" || !state.showObs) { $("#obs-note").textContent = ""; return; }
    const sp = spById(state.sp), e = evOf(sp, state.ev);
    const stage = e.stage || Object.entries(man.validate_map?.[sp.id] || {}).find(([, id]) => id === e.id)?.[0];
    let data;
    try { data = await (obsCache[sp.id] ||= fetchJSON(`obs/${sp.id}.json`)); } catch { $("#obs-note").textContent = ""; return; }
    const yr = yearNum();
    const recs = data.records.filter((r) => (!stage || r[4] === stage) && r[2] === yr);
    for (const r of recs) L.circleMarker([r[1], r[0]], { renderer: obsRenderer, radius: 3, weight: 1, color: "#111", fillColor: "#fff", fillOpacity: 0.9, interactive: false }).addTo(obsLayer);
    sizeObs();
    $("#obs-note").textContent = recs.length ? `${recs.length.toLocaleString()} iNaturalist ${stage ? stage.toLowerCase() + " " : ""}records from ${yr} shown as white dots.` : `No iNaturalist ${stage ? stage.toLowerCase() + " " : ""}records from ${yr} yet.`;
  }
  function sizeObs() {
    const r = Math.max(1.5, Math.min(5, 0.8 * 2 ** (map.getZoom() - 4)));
    obsLayer.eachLayer((m) => m.setRadius(r));
  }

  // ---------------------------------------------------------------- controls
  function fillSpecies() {
    const sel = $("#sp"), q = $("#sp-filter").value.trim().toLowerCase();
    sel.replaceChildren();
    for (const [grp, label] of [["insect", "Insects"], ["plant", "Plants"]]) {
      const og = h("optgroup", { label });
      for (const s of man.species.filter((x) => x.group === grp)) {
        if (q && !`${s.common} ${s.scientific}`.toLowerCase().includes(q)) continue;
        og.append(h("option", { value: s.id }, `${s.common} (${s.scientific})`));
      }
      if (og.children.length) sel.append(og);
    }
    if (![...sel.options].some((o) => o.value === state.sp) && sel.options.length) state.sp = sel.options[0].value;
    sel.value = state.sp;
  }
  function fillEvents() {
    const sp = spById(state.sp), sel = $("#ev");
    sel.replaceChildren(...sp.events.map((e) => h("option", { value: e.id }, e.label)));
    if (!evOf(sp, state.ev)) state.ev = sp.events[0].id;
    sel.value = state.ev;
  }
  function fillYears() {
    const sel = $("#year");
    sel.replaceChildren(h("option", { value: "current" }, `${curYear()} with forecast`));
    for (const y of Object.keys(man.years).sort().reverse()) sel.append(h("option", { value: y }, y));
    sel.value = state.year;
  }
  function syncControls() {
    for (const b of document.querySelectorAll("#mode button")) b.setAttribute("aria-pressed", String(b.dataset.mode === state.mode));
    $("#species-controls").hidden = state.mode !== "species";
    $("#today-controls").hidden = state.mode !== "today";
    for (const b of document.querySelectorAll("#views button")) b.setAttribute("aria-pressed", String(b.dataset.view === state.view));
    for (const b of document.querySelectorAll("#tviews button")) b.setAttribute("aria-pressed", String(b.dataset.tview === state.tview));
    for (const b of document.querySelectorAll("#units button")) b.setAttribute("aria-pressed", String(b.dataset.units === state.units));
    $("#driver").value = state.driver;
    $("#range-only").checked = state.rangeOnly;
    $("#show-obs").checked = state.showObs;
    $("#status-row").hidden = !(state.mode === "species" && state.view === "status");
    $("#stage-doy").value = state.stageDoy;
    $("#stage-src").value = state.stageSrc;
    $("#stage-date").textContent = state.stageSrc === "normal" ? M.fmtDoy(state.stageDoy) : M.fmtDoy(state.stageDoy, curYear());
    $("#ev").disabled = state.view === "status";
    $("#year").disabled = state.view === "status";
    renderModelInfo();
  }

  function renderModelInfo() {
    const box = $("#model-info");
    if (state.mode !== "species") { box.hidden = true; return; }
    box.hidden = false;
    const sp = spById(state.sp), e = evOf(sp, state.ev);
    const src = man.sources[e.source] || {};
    const val = validation[sp.id] || {};
    const parts = [h("p", { class: "label" }, "About this prediction")];
    const m = e.model;
    let how = "";
    if (m.type === "gdd") how = `${fmtDD(m.threshold)} ${ddUnit()} above ${state.units === "F" ? Math.round(m.base * 1.8 + 32) + " °F" : m.base + " °C"}${m.upper ? `, upper cutoff ${state.units === "F" ? Math.round(m.upper * 1.8 + 32) + " °F" : m.upper + " °C"}` : ""}, counted from ${typeof m.start === "string" ? fmtISO("2021-" + m.start) : `day length ${m.start.daylength} h`}.`;
    else if (m.type === "cold") how = `${Math.round(state.units === "F" ? m.threshold * 1.8 : m.threshold)} chilling ${ddUnit()} below ${state.units === "F" ? Math.round(m.base * 1.8 + 32) + " °F" : m.base + " °C"}, counted from ${typeof m.start === "string" ? fmtISO("2021-" + m.start) : `when days fall below ${m.start.daylength} h`}.`;
    else if (m.type === "precip") how = `At least ${(m.mm / 25.4).toFixed(1)} in. (${Math.round(m.mm)} mm) of rain in ${m.window} days, between ${fmtISO("2021-" + m.start)} and ${fmtISO("2021-" + m.end)}.`;
    else if (m.type === "photoperiod") how = `Day length reaches ${m.hours} h.`;
    else if (m.type === "shortday") how = `Days shorten to ${(16 - m.threshold).toFixed(1)} h of daylight (after midsummer).`;
    parts.push(h("p", {}, h("strong", {}, e.label), `: ${how}`));
    if (e.note) parts.push(h("p", { class: "note" }, e.note));
    if (e.stage && STAGE_ACTIONS[e.stage]) parts.push(h("p", { class: "action" }, STAGE_ACTIONS[e.stage]));
    const stage = e.stage || Object.keys(val).find((st) => val[st].kind === "literature" && val[st].event === e.id);
    const v = stage ? val[stage] : null;
    if (v && v.kind === "fitted" && v.test) {
      const t = v.test, t0 = v.test_null;
      const err = t.regions >= 5
        ? `In held-out years, predicted start and peak dates were off by ${t.region_error.toFixed(0)} days on average across ${t.regions} region-years, against ${t0.region_error.toFixed(0)} days for a single calendar date.`
        : `In held-out years the peak date was off by ${t.mae.toFixed(0)} days on average, against ${t0.mae.toFixed(0)} days for a single calendar date.`;
      parts.push(h("p", { class: "verdict" }, h("span", { class: "badge fitted" }, "Fitted"),
        ` from ${(v.n_train + v.n_test).toLocaleString()} iNaturalist ${stage.toLowerCase()} records (${v.driver_label}). ${err} ${Math.round(t.coverage * 100)}% of records fell between the predicted start and end.`));
    } else if (v && v.kind === "fitted") {
      parts.push(h("p", { class: "verdict" }, h("span", { class: "badge fitted" }, "Fitted"), ` from ${(v.n_train || 0).toLocaleString()} iNaturalist records; ${v.validated_by}.`));
    } else if (v && v.kind === "literature" && v.n) {
      parts.push(h("p", { class: "verdict" }, h("span", { class: "badge lit" }, "Published model"),
        ` Checked against ${v.n.toLocaleString()} iNaturalist ${stage.toLowerCase()} records: ${Math.round(v.early * 100)}% were photographed before the predicted date (${Math.round(v.early_7d * 100)}% more than a week before).`));
    } else {
      parts.push(h("p", { class: "verdict" }, h("span", { class: "badge lit" }, "Published model"), " Not yet checked against observations."));
    }
    parts.push(h("p", { class: "note" }, "Source: ", src.url ? h("a", { href: src.url }, src.cite) : src.cite || e.source));
    box.replaceChildren(...parts);
  }

  function wire() {
    for (const b of document.querySelectorAll("#mode button")) b.addEventListener("click", () => { state.mode = b.dataset.mode; refresh(); });
    $("#sp-filter").addEventListener("input", () => { fillSpecies(); fillEvents(); refresh(); });
    $("#sp").addEventListener("change", (e) => { state.sp = e.target.value; state.ev = null; fillEvents(); refresh(); });
    $("#ev").addEventListener("change", (e) => { state.ev = e.target.value; refresh(); });
    $("#year").addEventListener("change", (e) => { state.year = e.target.value; refresh(); });
    for (const b of document.querySelectorAll("#views button")) b.addEventListener("click", () => { state.view = b.dataset.view; refresh(); });
    for (const b of document.querySelectorAll("#tviews button")) b.addEventListener("click", () => { state.tview = b.dataset.tview; refresh(); });
    for (const b of document.querySelectorAll("#units button")) b.addEventListener("click", () => { state.units = b.dataset.units; refresh(); renderSide(); });
    $("#driver").addEventListener("change", (e) => { state.driver = e.target.value; refresh(); renderSide(); });
    $("#range-only").addEventListener("change", (e) => { state.rangeOnly = e.target.checked; refresh(); });
    $("#show-obs").addEventListener("change", (e) => { state.showObs = e.target.checked; drawObs(); writeHash(); });
    $("#stage-doy").addEventListener("input", (e) => { state.stageDoy = +e.target.value; syncControls(); throttled(); });
    $("#stage-src").addEventListener("change", (e) => { state.stageSrc = e.target.value; refresh(); });
    $("#play").addEventListener("click", togglePlay);
    $("#search").addEventListener("submit", onSearch);
    $("#locate").addEventListener("click", onLocate);
    $("#copy-link").addEventListener("click", async () => {
      const b = $("#copy-link");
      try { await navigator.clipboard.writeText(location.href); b.textContent = "Link copied"; } catch { b.textContent = "Copy the address bar"; }
      setTimeout(() => { b.textContent = "Copy link"; }, 2000);
    });
  }
  function refresh() { syncControls(); renderMap(); writeHash(); }
  let raf = null;
  const throttled = () => { if (!raf) raf = requestAnimationFrame(() => { raf = null; renderMap(); writeHash(); }); };
  let playTimer = null;
  function togglePlay() {
    const b = $("#play");
    if (playTimer) { clearInterval(playTimer); playTimer = null; b.textContent = "Play season"; return; }
    b.textContent = "Pause";
    if (state.stageDoy >= 340) state.stageDoy = 40;
    const reduce = matchMedia("(prefers-reduced-motion: reduce)").matches;
    playTimer = setInterval(() => {
      state.stageDoy += reduce ? 14 : 3;
      if (state.stageDoy > 350) { state.stageDoy = 350; togglePlay(); }
      syncControls(); throttled();
    }, reduce ? 900 : 120);
  }

  // ---------------------------------------------------------------- location panel
  async function selectLocation(lat, lon, name = null) {
    state.loc = { lat, lon, name };
    marker.setLatLng([lat, lon]);
    writeHash();
    renderSide();
  }

  async function pointData(lat, lon) {
    const c = cellOf(lat, lon);
    if (!c) return null;
    const T = man.tile, name = `${Math.floor(c.i / T)}_${Math.floor(c.j / T)}`;
    if (!man.tiles.includes(name)) return null;
    const a = await fetchBin(`tiles/${name}.bin`);
    const cell = (c.i % T) * T + (c.j % T), S = T * T, nE = man.tile_events.length;
    const events = {};
    man.tile_events.forEach((key, k) => {
      const { doy, src } = M.decodeCurrent(a[(k * 3 + 2) * S + cell]);
      events[key] = { mean: a[(k * 3) * S + cell] / 10, sd: a[(k * 3 + 1) * S + cell] / 10, doy, src };
    });
    const range = {};
    man.tile_species.forEach((id, k) => { range[id] = a[(nE * 3 + k) * S + cell]; });
    const general = {};
    for (const d of man.drivers) {
      const g = await generalLayer(d.id), n = N();
      general[d.id] = g[c.idx] === 65535 ? null : { value: g[c.idx], normal: g[n + c.idx], ahead: g[2 * n + c.idx] === 65535 ? null : g[2 * n + c.idx] - 1000 };
    }
    return { events, range, general, idx: c.idx };
  }

  let sideToken = null;
  async function renderSide() {
    const token = (sideToken = {});
    const { lat, lon, name } = state.loc;
    $("#loc-name").textContent = name || "Selected location";
    $("#loc-coords").textContent = `${Math.abs(lat).toFixed(2)}° N, ${Math.abs(lon).toFixed(2)}° W`;
    let p = null;
    try { p = await pointData(lat, lon); } catch { p = null; }
    if (token !== sideToken) return;
    const body = $("#side-body");
    if (!p) { body.replaceChildren(h("p", { class: "note" }, "This tool covers the contiguous United States. Pick a point inside it.")); return; }
    const td = todayDoy(), y = curYear();
    // Degree-days here
    const g = p.general[state.driver];
    const ddBox = h("div", { class: "dd-now" });
    if (g) {
      const ahead = g.ahead;
      ddBox.append(
        h("div", { class: "stat" }, h("span", { class: "label" }, `${man.drivers.find((d) => d.id === state.driver).label}, to date`), h("strong", {}, `${fmtDD(g.value)}`), h("span", { class: "unit" }, ddUnit())),
        h("div", { class: "stat" }, h("span", { class: "label" }, "Normal to date"), h("strong", {}, fmtDD(g.normal)), h("span", { class: "unit" }, ddUnit())),
        h("div", { class: "stat" }, h("span", { class: "label" }, "Season"), h("strong", { class: ahead > 0 ? "ahead" : ahead < 0 ? "behind" : "" }, ahead == null ? "–" : ahead === 2000 ? "Past" : ahead === 0 ? "On time" : `${Math.abs(ahead)} d`),
          h("span", { class: "unit" }, ahead === 2000 ? "a normal full year's total" : ahead > 0 ? "ahead of normal" : ahead < 0 ? "behind normal" : "")));
    }
    // Upcoming
    const near = (id) => !state.rangeOnly || p.range[id] !== 0; // 2 = no range data yet
    const rows = [];
    for (const sp of man.species) {
      if (!near(sp.id)) continue;
      for (const e of sp.events) {
        const x = p.events[`${sp.id}.${e.id}`];
        if (!x || !x.doy) continue;
        const dd = x.doy - td;
        if (dd >= -7 && dd <= 21) rows.push({ sp, e, x, dd });
      }
    }
    rows.sort((a, b) => a.x.doy - b.x.doy);
    const upcoming = h("ol", { class: "upcoming" }, rows.length ? rows.map(({ sp, e, x, dd }) =>
      h("li", { class: dd < 0 ? "past" : "" },
        h("span", { class: "when" }, M.fmtDoy(x.doy, y), h("small", {}, M.relDays(dd))),
        h("button", { type: "button", class: "what", onclick: () => { state.mode = "species"; state.sp = sp.id; state.ev = e.id; state.year = "current"; state.view = "date"; $("#sp-filter").value = ""; fillSpecies(); fillEvents(); refresh(); } },
          h("strong", {}, sp.common), ` ${e.label.toLowerCase()}`),
        h("span", { class: `src src-${x.src}` }, srcWord(x.src)))) : [h("li", { class: "empty" }, "Nothing modeled within a week either side of today here.")]);
    // Season calendar
    const cal = calendar(p, near, td, y);
    body.replaceChildren(
      ddBox,
      h("h3", { class: "label" }, "Coming up here (last 7 days to next 3 weeks)"), upcoming,
      h("h3", { class: "label" }, `Season at a glance, ${y}`), cal,
      h("p", { class: "note" }, state.rangeOnly ? "Lists species with iNaturalist records within about 25 km. Untick “Only where recorded” to see every modeled species." : "All modeled species, recorded nearby or not."));
  }

  function calendar(p, near, td, y) {
    const ns = "http://www.w3.org/2000/svg";
    const s = (t, a, txt) => { const e = document.createElementNS(ns, t); for (const [k, v] of Object.entries(a)) e.setAttribute(k, v); if (txt != null) e.textContent = txt; return e; };
    const list = man.species.filter((sp) => near(sp.id) && sp.events.some((e) => p.events[`${sp.id}.${e.id}`].mean > 0 || p.events[`${sp.id}.${e.id}`].doy > 0));
    const W = 380, L = 112, R = 8, rowH = 20, top = 18, H = top + list.length * rowH + 6;
    const X = (d) => L + ((d - 1) / 365) * (W - L - R);
    const svg = s("svg", { viewBox: `0 0 ${W} ${H}`, role: "img", "aria-label": "Season calendar for species at this location" });
    for (let m = 0; m < 12; m++) {
      const d = M.dateToDoy(new Date(Date.UTC(2021, m, 1)));
      svg.append(s("line", { x1: X(d), x2: X(d), y1: top - 4, y2: H, class: "gridline" }));
      if (m % 2 === 0) svg.append(s("text", { x: X(d) + 2, y: 10, class: "tick" }, M.MONTHS[m]));
    }
    svg.append(s("line", { x1: X(td), x2: X(td), y1: top - 6, y2: H, class: "today-line" }));
    list.forEach((sp, r) => {
      const yy = top + r * rowH + rowH / 2;
      const name = s("text", { x: 0, y: yy + 4, class: "row-label" }, sp.common.length > 19 ? sp.common.slice(0, 18) + "…" : sp.common);
      svg.append(name);
      const evs = sp.events.map((e) => ({ e, x: p.events[`${sp.id}.${e.id}`] }));
      for (const { e, x } of evs) {
        if (x.mean > 0) {
          if (x.sd > 0) svg.append(s("line", { x1: X(x.mean - x.sd), x2: X(x.mean + x.sd), y1: yy, y2: yy, class: "sd-bar" }));
          svg.append(s("line", { x1: X(x.mean), x2: X(x.mean), y1: yy - 5, y2: yy + 5, class: "normal-tick" }));
        }
        if (x.doy > 0) {
          const c = s("circle", { cx: X(x.doy), cy: yy, r: 3.6, class: `dot src-${x.src}` });
          c.append(s("title", {}, `${sp.common}: ${e.label} ${M.fmtDoy(x.doy, y)} (${srcWord(x.src)})${x.mean ? `; normal ${M.fmtDoy(x.mean)}` : ""}`));
          svg.append(c);
        }
      }
    });
    if (!list.length) return h("p", { class: "note" }, "No modeled species recorded near this location.");
    return h("div", { class: "cal" }, svg, h("ul", { class: "key" },
      h("li", {}, h("i", { class: "k-normal" }), "Normal ± 1 SD"), h("li", {}, h("i", { class: "k-dot src-1" }), "Observed"),
      h("li", {}, h("i", { class: "k-dot src-3" }), "Forecast"), h("li", {}, h("i", { class: "k-dot src-4" }), "Projected")));
  }

  // ---------------------------------------------------------------- search / locate / hash
  async function onSearch(e) {
    e.preventDefault();
    const text = $("#q").value.trim(), msg = $("#search-msg");
    msg.textContent = "";
    const ll = text.match(/^\s*(-?\d+(?:\.\d+)?)\s*[, ]\s*(-?\d+(?:\.\d+)?)\s*$/);
    if (ll) { selectLocation(+ll[1], +ll[2], `${(+ll[1]).toFixed(2)}, ${(+ll[2]).toFixed(2)}`); map.setView([+ll[1], +ll[2]], Math.max(6, map.getZoom())); return; }
    if (!text) return;
    try {
      const r = await fetch(`https://nominatim.openstreetmap.org/search?format=jsonv2&limit=1&countrycodes=us&q=${encodeURIComponent(text)}`);
      const j = await r.json();
      if (!j.length) { msg.textContent = `No place found for “${text}”. Try a town and state, or “lat, lon”.`; return; }
      selectLocation(+j[0].lat, +j[0].lon, j[0].display_name.split(",").slice(0, 2).join(","));
      map.setView([+j[0].lat, +j[0].lon], Math.max(6, map.getZoom()));
    } catch { msg.textContent = "Place search is not reachable here. Enter coordinates as “lat, lon”, or click the map."; }
  }
  function onLocate() {
    const msg = $("#search-msg");
    if (!navigator.geolocation) { msg.textContent = "This browser does not share location."; return; }
    navigator.geolocation.getCurrentPosition((p) => { selectLocation(p.coords.latitude, p.coords.longitude, "Your location"); map.setView([p.coords.latitude, p.coords.longitude], 7); },
      () => { msg.textContent = "Location access was not allowed. Click the map or search instead."; }, { timeout: 10000 });
  }
  function writeHash() {
    const p = new URLSearchParams({ m: state.mode, u: state.units, at: `${state.loc.lat.toFixed(3)},${state.loc.lon.toFixed(3)}` });
    if (state.mode === "today") { p.set("d", state.driver); p.set("v", state.tview); }
    else { p.set("sp", state.sp); p.set("e", state.ev); p.set("y", state.year); p.set("v", state.view); if (state.view === "status") { p.set("doy", state.stageDoy); p.set("s", state.stageSrc); } }
    if (!state.rangeOnly) p.set("r", "0");
    try { history.replaceState(null, "", "#" + p.toString()); } catch { /* sandboxed */ }
  }
  function readHash() {
    if (!location.hash.includes("=")) return;
    const p = new URLSearchParams(location.hash.slice(1));
    if (["today", "species"].includes(p.get("m"))) state.mode = p.get("m");
    if (["F", "C"].includes(p.get("u"))) state.units = p.get("u");
    if (p.get("d") && man.drivers.some((d) => d.id === p.get("d"))) state.driver = p.get("d");
    if (p.get("sp") && spById(p.get("sp"))) state.sp = p.get("sp");
    if (p.get("e")) state.ev = p.get("e");
    if (p.get("y") && (p.get("y") === "current" || man.years[p.get("y")])) state.year = p.get("y");
    if (p.get("v")) { if (state.mode === "today" && TODAY_VIEWS[p.get("v")]) state.tview = p.get("v"); if (state.mode === "species" && SPECIES_VIEWS[p.get("v")]) state.view = p.get("v"); }
    if (p.get("doy")) state.stageDoy = Math.min(366, Math.max(1, +p.get("doy")));
    if (["current", "normal"].includes(p.get("s"))) state.stageSrc = p.get("s");
    if (p.get("r") === "0") state.rangeOnly = false;
    const at = (p.get("at") || "").split(",").map(Number);
    if (at.length === 2 && at.every(Number.isFinite)) state.loc = { lat: at[0], lon: at[1], name: null };
  }

  // ---------------------------------------------------------------- boot
  async function boot() {
    try { man = await fetchJSON("manifest.json"); } catch {
      $("#map-title").textContent = "The map data could not be loaded. Reload the page to try again."; return;
    }
    try { validation = await fetchJSON("validation.json"); } catch { validation = {}; }
    man.validate_map = {};
    for (const [sp, st] of Object.entries(validation)) for (const [stage, v] of Object.entries(st)) if (v.kind === "literature") (man.validate_map[sp] ||= {})[stage] = v.event;
    state.stageDoy = todayDoy();
    readHash();
    if (!spById(state.sp)) state.sp = man.species[0].id;
    const c = cur();
    $("#status").replaceChildren(
      h("span", {}, h("b", {}, M.fmtDoy(c.today_doy, c.year)), ` ${c.year}`),
      h("span", {}, `Observed through ${fmtISO(c.obs_through)}`),
      h("span", {}, c.nbm_cycle ? `NWS forecast to ${M.fmtDoy(c.forecast_through_doy, c.year)}` : "Forecast unavailable today"),
      h("span", {}, "Normals after that"));
    $("#driver").replaceChildren(...man.drivers.map((d) => h("option", { value: d.id }, d.label)));
    fillSpecies(); fillEvents(); fillYears();
    if (/claude/.test(location.hostname)) $("#copy-link").hidden = true;
    try {
      const b = await fetchJSON("boundaries.json");
      boundaries = b.map((f) => { const o = []; for (let k = 0; k < f.length; k += 2) o.push([f[k + 1], f[k]]); return o; });
    } catch { /* optional */ }
    const g = await generalLayer(man.drivers[0].id), n = N();
    landMask = new Uint8Array(n);
    for (let i = 0; i < n; i++) landMask[i] = g[i] !== 65535 ? 1 : 0;
    initMap();
    wire();
    refresh();
    renderSide();
    matchMedia("(prefers-color-scheme: dark)").addEventListener("change", restyle);
    new MutationObserver(restyle).observe(document.documentElement, { attributes: true, attributeFilter: ["data-theme"] });
  }
  boot();
})();
