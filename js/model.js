/* Date helpers and data decoding shared by the app (and the Node tests). */
(function (root) {
  "use strict";
  const MONTHS = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"];
  const isLeap = (y) => (y % 4 === 0 && y % 100 !== 0) || y % 400 === 0;
  const doyToDate = (doy, y) => new Date(Date.UTC(y, 0, Math.round(doy)));
  const dateToDoy = (d) => Math.round((Date.UTC(d.getUTCFullYear(), d.getUTCMonth(), d.getUTCDate()) - Date.UTC(d.getUTCFullYear(), 0, 0)) / 864e5);
  /** "Apr 23"; normals use a non-leap reference year. */
  function fmtDoy(doy, y = 2021) {
    const d = doyToDate(doy, y);
    return `${MONTHS[d.getUTCMonth()]} ${d.getUTCDate()}`;
  }
  /** Current-year layer value -> {doy, src}; src 1 observed, 2 provisional, 3 forecast, 4 projected. */
  const decodeCurrent = (v) => ({ doy: v & 511, src: v >> 9 });
  const SOURCES = { 1: "observed", 2: "provisional", 3: "forecast", 4: "projected" };
  /** Signed day difference phrased for people. */
  function relDays(n) {
    if (n === 0) return "today";
    if (n === 1) return "tomorrow";
    if (n === -1) return "yesterday";
    return n > 0 ? `in ${n} days` : `${-n} days ago`;
  }
  const api = { MONTHS, isLeap, doyToDate, dateToDoy, fmtDoy, decodeCurrent, SOURCES, relDays };
  if (typeof module !== "undefined" && module.exports) module.exports = api;
  else root.PhenoModel = api;
})(this);
