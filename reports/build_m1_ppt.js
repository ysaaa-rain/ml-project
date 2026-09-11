const pptxgen = require('C:/Users/36046/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/pptxgenjs');
const path = require('path');

const pptx = new pptxgen();
pptx.layout = 'LAYOUT_WIDE';
pptx.author = 'MLG04 project team';
pptx.subject = '启动子调控模式发现：第一阶段问题建模与基础流程验证';
pptx.title = '启动子调控模式发现｜第一阶段进度汇报';
pptx.company = 'MLG04';
pptx.lang = 'zh-CN';
pptx.theme = {
  headFontFace: 'Microsoft YaHei',
  bodyFontFace: 'Microsoft YaHei',
  lang: 'zh-CN',
};
pptx.defineSlideMaster({
  title: 'MASTER',
  background: { color: 'F6FAFC' },
  objects: [
    { line: { x: 0.55, y: 7.12, w: 12.2, h: 0, line: { color: 'D8E5EA', width: 0.8 } } },
  ],
  slideNumber: { x: 12.3, y: 7.10, color: '6B8190', fontFace: 'Aptos', fontSize: 8 },
});

const W = 13.333;
const H = 7.5;
const C = {
  navy: '173650',
  ink: '24455D',
  muted: '6B8190',
  teal: '2AA198',
  tealDark: '16766F',
  mint: 'E1F3EE',
  blue: 'E8F2F8',
  orange: 'FCEAD8',
  orangeDark: 'C96E27',
  red: 'FBE7E5',
  redDark: 'C74A45',
  gray: 'EEF3F5',
  line: 'D4E1E7',
  white: 'FFFFFF',
};

const ROOT = path.resolve(__dirname, '..');
const IMG_CLEAN = path.join(ROOT, 'reports', 'M1证据_数据清洗与分组.png');
const IMG_B0 = path.join(ROOT, 'reports', 'M1证据_B0已知元件扫描.png');
const IMG_ENV = path.join(ROOT, 'reports', 'M1证据_可复现环境检查.png');
const OUT = path.join(ROOT, 'reports', 'MLG04_第一阶段进度汇报_问题建模与基础流程验证.pptx');

function tx(slide, text, x, y, w, h, opts = {}) {
  slide.addText(text, {
    x, y, w, h,
    fontFace: opts.fontFace || 'Microsoft YaHei',
    fontSize: opts.fontSize || 16,
    color: opts.color || C.ink,
    bold: opts.bold || false,
    margin: opts.margin === undefined ? 0 : opts.margin,
    breakLine: false,
    fit: 'shrink',
    valign: opts.valign || 'mid',
    align: opts.align || 'left',
    paraSpaceAfterPt: opts.paraSpaceAfterPt || 0,
    italic: opts.italic || false,
    bullet: opts.bullet,
    transparency: opts.transparency,
    charSpacing: opts.charSpacing,
  });
}

function rect(slide, x, y, w, h, fill, radius = 0.12, line = fill) {
  slide.addShape(pptx.ShapeType.roundRect, {
    x, y, w, h,
    rectRadius: radius,
    fill: { color: fill },
    line: { color: line, width: 0.8 },
  });
}

function line(slide, x, y, w, h, color = C.line, width = 1.2, dash = 'solid') {
  slide.addShape(pptx.ShapeType.line, { x, y, w, h, line: { color, width, dashType: dash, beginArrowType: 'none', endArrowType: 'none' } });
}

function top(slide, kicker, title, subtitle = '') {
  tx(slide, kicker.toUpperCase(), 0.62, 0.28, 2.4, 0.24, { fontFace: 'Aptos', fontSize: 9, bold: true, color: C.tealDark, charSpacing: 1.2 });
  tx(slide, title, 0.62, 0.58, 11.7, 0.48, { fontSize: 25, bold: true, color: C.navy });
  if (subtitle) tx(slide, subtitle, 0.64, 1.08, 11.9, 0.28, { fontSize: 11.5, color: C.muted });
}

function pill(slide, label, x, y, w, fill = C.mint, color = C.tealDark) {
  rect(slide, x, y, w, 0.31, fill, 0.15, fill);
  tx(slide, label, x, y + 0.005, w, 0.28, { fontSize: 9.5, bold: true, color, align: 'center' });
}

function bulletList(slide, items, x, y, w, lineH = 0.38, fontSize = 14, color = C.ink, bulletColor = C.teal) {
  items.forEach((item, i) => {
    const yy = y + i * lineH;
    slide.addShape(pptx.ShapeType.ellipse, { x, y: yy + (lineH - 0.09) / 2, w: 0.09, h: 0.09, fill: { color: bulletColor }, line: { color: bulletColor } });
    tx(slide, item, x + 0.20, yy, w - 0.20, lineH, { fontSize, color });
  });
}

function metric(slide, value, label, x, y, w, fill = C.blue, valueColor = C.navy) {
  rect(slide, x, y, w, 0.95, fill, 0.14, fill);
  tx(slide, value, x + 0.16, y + 0.10, w - 0.32, 0.40, { fontFace: 'Aptos Display', fontSize: 24, bold: true, color: valueColor });
  tx(slide, label, x + 0.16, y + 0.58, w - 0.32, 0.20, { fontSize: 9.5, color: C.muted });
}

function arrow(slide, x, y, w, color = C.teal) {
  slide.addShape(pptx.ShapeType.chevron, { x, y, w, h: 0.25, fill: { color }, line: { color } });
}

// 1. Title
{
  const s = pptx.addSlide('MASTER');
  s.background = { color: C.navy };
  s.addShape(pptx.ShapeType.rect, { x: 0, y: 0, w: W, h: H, fill: { color: C.navy }, line: { color: C.navy } });
  tx(s, 'MLG04 / M1', 0.72, 0.55, 2.1, 0.26, { fontFace: 'Aptos', fontSize: 11, bold: true, color: '8ED6CB', charSpacing: 1.5 });
  tx(s, '启动子调控模式发现', 0.72, 1.55, 6.8, 0.72, { fontSize: 34, bold: true, color: C.white });
  tx(s, '第一阶段进度汇报', 0.75, 2.48, 5.4, 0.48, { fontSize: 23, bold: true, color: 'DFF5F1' });
  tx(s, '问题建模与基础流程验证', 0.76, 3.06, 5.7, 0.34, { fontSize: 15, color: 'A9C5D0' });
  line(s, 0.76, 3.75, 2.3, 0, C.teal, 3);
  tx(s, '课程项目：机器学习综合实践', 0.76, 6.55, 4.4, 0.26, { fontSize: 10.5, color: 'A9C5D0' });
  // Minimal motif illustration
  tx(s, 'DNA sequence', 7.42, 1.10, 1.6, 0.25, { fontFace: 'Aptos', fontSize: 10, color: 'A9C5D0' });
  const bases = ['A', 'C', 'G', 'T', 'T', 'G', 'A', 'C', 'A', 'T', 'A', 'A', 'T'];
  bases.forEach((b, i) => {
    const xx = 7.4 + i * 0.37;
    const active = i >= 5 && i <= 10;
    s.addShape(pptx.ShapeType.roundRect, { x: xx, y: 2.18, w: 0.29, h: 0.42, rectRadius: 0.06, fill: { color: active ? C.teal : '31536B' }, line: { color: active ? C.teal : '31536B' } });
    tx(s, b, xx, 2.20, 0.29, 0.34, { fontFace: 'Aptos', fontSize: 15, bold: true, color: C.white, align: 'center' });
  });
  line(s, 7.52, 2.95, 2.0, 0, '6FAEA8', 1.5);
  tx(s, 'motif', 8.22, 3.08, 1.0, 0.25, { fontFace: 'Aptos', fontSize: 11, bold: true, color: '8ED6CB', align: 'center' });
  line(s, 10.15, 2.95, 1.65, 0, '6FAEA8', 1.5);
  tx(s, '位置 / 间距 / 组合', 8.06, 4.05, 3.3, 0.32, { fontSize: 14, color: 'DFF5F1', align: 'center' });
  tx(s, '从序列中找规律，再验证规律是否可靠', 7.25, 5.25, 4.7, 0.34, { fontSize: 14, color: 'A9C5D0', align: 'center' });
}

// 2. Problem
{
  const s = pptx.addSlide('MASTER');
  top(s, '01 / 问题', '项目要解决什么问题', '目标：从启动子序列中发现调控模式');
  rect(s, 0.62, 1.63, 3.05, 4.78, C.navy, 0.18, C.navy);
  tx(s, '启动子', 0.92, 1.98, 2.4, 0.48, { fontSize: 28, bold: true, color: C.white });
  tx(s, '可以理解为基因前面的\n“开关区域”', 0.92, 2.63, 2.35, 0.80, { fontSize: 19, bold: true, color: 'DFF5F1' });
  line(s, 0.94, 3.77, 1.25, 0, C.teal, 3);
  tx(s, '我们关注序列中的\n重复片段与排列规律', 0.92, 4.05, 2.35, 0.76, { fontSize: 15, color: 'A9C5D0' });
  tx(s, '无监督模式发现', 0.92, 5.70, 2.1, 0.28, { fontSize: 13, bold: true, color: '8ED6CB' });

  const cards = [
    ['01', '哪些短序列反复出现？', 'motif'],
    ['02', '它们出现在什么位置？', '位置'],
    ['03', '是否存在固定间距或组合？', '排列'],
    ['04', '换组别后还能复现吗？', '验证'],
  ];
  cards.forEach((c, i) => {
    const x = 4.08 + (i % 2) * 4.15;
    const y = 1.66 + Math.floor(i / 2) * 1.72;
    rect(s, x, y, 3.75, 1.36, i === 0 ? C.mint : C.white, 0.16, C.line);
    tx(s, c[0], x + 0.25, y + 0.22, 0.45, 0.28, { fontFace: 'Aptos', fontSize: 12, bold: true, color: C.tealDark });
    tx(s, c[1], x + 0.25, y + 0.54, 3.15, 0.36, { fontSize: 15, bold: true, color: C.navy });
    tx(s, c[2], x + 0.25, y + 1.00, 1.2, 0.18, { fontFace: 'Aptos', fontSize: 10, color: C.muted });
  });
  rect(s, 4.08, 5.30, 8.00, 0.78, C.orange, 0.15, C.orange);
  tx(s, '本项目不是：', 4.35, 5.49, 1.25, 0.25, { fontSize: 12, bold: true, color: C.orangeDark });
  tx(s, '启动子分类  ·  强度预测  ·  序列生成', 5.62, 5.49, 5.95, 0.25, { fontSize: 14, bold: true, color: C.ink });
}

// 3. Research questions
{
  const s = pptx.addSlide('MASTER');
  top(s, '02 / 建模', '第一阶段核心：把问题变成可计算问题', '一个开放的生物学问题，被拆成三个可以统计和验证的问题');
  const qs = [
    ['RQ1', '发现的 motif 是否对应已知调控元件？', '比较 -10 box、-35 box、UP element 等已知元件。', C.mint, C.tealDark],
    ['RQ2', '不同 σ 因子启动子是否有不同规律？', '比较 motif 频率、位置、间距和组合。', C.blue, C.navy],
    ['RQ3', '规律能否跨物种、跨数据集复现？', '区分稳定规律、物种特异规律、数据偏差和随机假象。', C.orange, C.orangeDark],
  ];
  qs.forEach((q, i) => {
    const y = 1.70 + i * 1.48;
    rect(s, 0.82, y, 11.70, 1.18, q[3], 0.18, q[3]);
    tx(s, q[0], 1.10, y + 0.25, 0.80, 0.30, { fontFace: 'Aptos', fontSize: 17, bold: true, color: q[4] });
    tx(s, q[1], 2.05, y + 0.20, 5.95, 0.34, { fontSize: 18, bold: true, color: C.navy });
    tx(s, q[2], 2.05, y + 0.66, 9.50, 0.25, { fontSize: 12.5, color: C.ink });
  });
  rect(s, 0.82, 6.22, 11.70, 0.48, C.navy, 0.13, C.navy);
  tx(s, '第一阶段最重要的工作：把“想研究什么”拆成后续程序能够计算、统计和验证的对象。', 1.10, 6.31, 11.12, 0.24, { fontSize: 12.5, color: C.white, align: 'center' });
}

// 4. Data modeling and evidence
{
  const s = pptx.addSlide('MASTER');
  top(s, '03 / 数据', '数据如何被建模', '先把每条启动子统一成标准记录，再进入 motif 分析');
  rect(s, 0.62, 1.55, 4.35, 4.92, C.white, 0.17, C.line);
  tx(s, '一条标准记录包含什么？', 0.94, 1.86, 3.5, 0.30, { fontSize: 17, bold: true, color: C.navy });
  pill(s, '核心字段', 0.94, 2.32, 1.05, C.mint, C.tealDark);
  tx(s, '`sequence_id`   `species`\n`source_dataset`   `sequence`', 0.94, 2.75, 3.55, 0.66, { fontFace: 'Aptos', fontSize: 15, bold: true, color: C.ink });
  pill(s, '可选字段', 0.94, 3.68, 1.05, C.blue, C.navy);
  tx(s, '`sigma_factor_type`\n`evidence_level`   `promoter_strength`', 0.94, 4.10, 3.55, 0.75, { fontFace: 'Aptos', fontSize: 14, color: C.ink });
  line(s, 0.94, 5.18, 3.55, 0, C.line, 1);
  tx(s, '处理规则：统一格式 · 过滤问题记录\n去重 · discovery / validation 分组\n生成随机背景作为对照', 0.94, 5.42, 3.55, 0.64, { fontSize: 12.5, color: C.muted });
  // Screenshot panel
  rect(s, 5.25, 1.55, 7.45, 4.92, C.white, 0.17, C.line);
  tx(s, '真实运行截图：数据清洗与分组', 5.55, 1.83, 3.8, 0.25, { fontSize: 15, bold: true, color: C.navy });
  s.addImage({ path: IMG_CLEAN, x: 5.52, y: 2.18, w: 6.90, h: 3.88 });
  tx(s, 'demo：8 条输入 → 6 条有效；1 条高 N，1 条重复', 5.55, 6.16, 6.65, 0.22, { fontSize: 10.5, color: C.redDark, italic: true });
}

// 5. Model framework
{
  const s = pptx.addSlide('MASTER');
  top(s, '04 / 模型', '已建立的模型框架', '从已知元件基线开始，逐步扩展到 de novo 发现和独立验证');
  const models = [
    ['B0', '已知元件\nPWM 扫描', '已实现', C.mint, C.tealDark],
    ['B1', 'MEME\nde novo', '已设计', C.blue, C.navy],
    ['B2', 'STREME /\nDREME', '已设计', C.blue, C.navy],
    ['B3', '按 sigma 因子\n分组发现', '已设计', C.orange, C.orangeDark],
    ['B4', '跨物种 /\n跨数据源', '已设计', C.orange, C.orangeDark],
  ];
  models.forEach((m, i) => {
    const x = 0.75 + i * 2.42;
    rect(s, x, 1.78, 2.05, 2.32, m[3], 0.18, m[3]);
    tx(s, m[0], x + 0.20, 2.02, 0.55, 0.30, { fontFace: 'Aptos', fontSize: 16, bold: true, color: m[4] });
    tx(s, m[1], x + 0.20, 2.54, 1.62, 0.72, { fontSize: 17, bold: true, color: C.navy, valign: 'top' });
    pill(s, m[2], x + 0.20, 3.54, 0.88, m[2] === '已实现' ? C.teal : C.white, m[2] === '已实现' ? C.white : m[4]);
    if (i < models.length - 1) arrow(s, x + 2.13, 2.84, 0.25, i === 0 ? C.teal : C.line);
  });
  rect(s, 0.75, 4.65, 11.78, 1.02, C.white, 0.16, C.line);
  tx(s, '已知元件', 1.05, 4.94, 1.10, 0.24, { fontSize: 13, bold: true, color: C.navy });
  tx(s, '`TTGACA`  =  -35 box', 2.32, 4.91, 2.42, 0.28, { fontFace: 'Aptos', fontSize: 13, color: C.ink });
  tx(s, '`TATAAT`  =  -10 box', 4.82, 4.91, 2.42, 0.28, { fontFace: 'Aptos', fontSize: 13, color: C.ink });
  line(s, 7.48, 4.83, 0, 0.62, C.line, 1);
  tx(s, 'N0 随机背景：已实现', 7.78, 4.91, 2.10, 0.28, { fontSize: 13, color: C.tealDark, bold: true });
  tx(s, 'N1 二核苷酸背景：待实现', 10.00, 4.91, 2.15, 0.28, { fontSize: 12, color: C.muted });
  rect(s, 0.75, 6.05, 11.78, 0.56, C.navy, 0.12, C.navy);
  tx(s, '当前主线：B0 基线 → MEME / STREME 发现 → 位置与显著性 → 分组和跨数据集验证', 1.02, 6.19, 11.22, 0.22, { fontSize: 12.5, color: C.white, align: 'center' });
}

// 6. Demo results
{
  const s = pptx.addSlide('MASTER');
  top(s, '05 / 结果', 'Demo 跑通结果', '8 条 demo 数据跑通了完整基础流程；结果用于验证代码链路，不是正式生物学结论');
  rect(s, 0.62, 1.55, 3.42, 4.98, C.white, 0.17, C.line);
  tx(s, '清洗结果', 0.94, 1.86, 2.4, 0.30, { fontSize: 17, bold: true, color: C.navy });
  metric(s, '8 → 6', '输入记录 → 有效记录', 0.94, 2.35, 2.78, C.mint, C.tealDark);
  metric(s, '1 + 1', '高 N 含量 + 重复序列', 0.94, 3.48, 2.78, C.orange, C.orangeDark);
  metric(s, '48 / 12', '-10 box / -35 box 命中次数', 0.94, 4.61, 2.78, C.blue, C.navy);
  tx(s, '7/7 自动化测试通过', 0.94, 5.92, 2.78, 0.25, { fontSize: 13, bold: true, color: C.tealDark, align: 'center' });
  // B0 screenshot
  rect(s, 4.35, 1.55, 8.35, 4.98, C.white, 0.17, C.line);
  tx(s, '真实运行截图：B0 已知元件扫描', 4.67, 1.83, 4.1, 0.25, { fontSize: 15, bold: true, color: C.navy });
  s.addImage({ path: IMG_B0, x: 4.62, y: 2.18, w: 7.80, h: 4.39 });
  tx(s, '说明：截图中的 p_value_proxy 只是排序代理，不是正式 FIMO 显著性。', 4.68, 6.68, 7.25, 0.20, { fontSize: 9.5, color: C.redDark, italic: true });
}

// 7. Conclusion and next step
{
  const s = pptx.addSlide('MASTER');
  top(s, '06 / 结论', '当前结论与下一阶段', '第一阶段完成：问题建模 + 基础流程验证');
  rect(s, 0.62, 1.55, 5.35, 4.98, C.mint, 0.17, C.mint);
  tx(s, '已经完成', 0.96, 1.87, 2.2, 0.32, { fontSize: 19, bold: true, color: C.tealDark });
  bulletList(s, [
    '研究问题定义（RQ1–RQ3）',
    '数据结构、清洗规则与分组',
    'B0 已知元件 PWM 基线',
    '随机背景对照与 EDA 框架',
    '输入哈希、环境和结果记录',
    '小样本 demo 验证，测试 7/7 通过',
  ], 0.98, 2.42, 4.45, 0.48, 13, C.ink, C.teal);
  rect(s, 6.25, 1.55, 6.45, 2.25, C.orange, 0.17, C.orange);
  tx(s, '下一阶段：进入真实数据实验', 6.58, 1.87, 5.4, 0.30, { fontSize: 18, bold: true, color: C.orangeDark });
  const steps = ['数据审计', '确定窗口', '安装 MEME Suite', '运行 MEME / FIMO', '比较与验证'];
  steps.forEach((st, i) => {
    const x = 6.58 + i * 1.14;
    rect(s, x, 2.58, 0.98, 0.58, C.white, 0.12, C.white);
    tx(s, st, x + 0.03, 2.72, 0.92, 0.20, { fontSize: 9.5, bold: true, color: C.ink, align: 'center' });
    if (i < steps.length - 1) arrow(s, x + 1.00, 2.75, 0.12, C.orangeDark);
  });
  tx(s, '当前正式 motif 工具尚未安装，不能宣称真实实验已经完成。', 6.58, 3.34, 5.45, 0.25, { fontSize: 11.5, color: C.redDark, italic: true });
  rect(s, 6.25, 4.10, 6.45, 2.43, C.white, 0.17, C.line);
  tx(s, '环境检查：Python 管线可运行，MEME Suite 待补齐', 6.58, 4.34, 5.5, 0.25, { fontSize: 14, bold: true, color: C.navy });
  tx(s, 'Python 3.13.5\nseed 20260911\n5 个主线工具当前未检测到', 6.58, 4.88, 2.05, 0.92, { fontSize: 11.5, color: C.ink, valign: 'top' });
  // Keep the evidence screenshot at its original 16:9 aspect ratio.
  s.addImage({ path: IMG_ENV, x: 8.82, y: 4.66, w: 3.52, h: 1.98 });
  tx(s, '阶段结论：完成的是“如何研究这个问题”的建模工作；下一阶段才回答真实数据中的 motif 问题。', 0.96, 6.62, 11.5, 0.24, { fontSize: 12.5, bold: true, color: C.navy, align: 'center' });
}

pptx.writeFile({ fileName: OUT });
console.log(OUT);
