# NanoRNAMod Benchmarking Workflow - Figure Description for AI Image Generation

## Overview
A scientific workflow diagram for RNA modification detection benchmarking from Oxford Nanopore Direct RNA Sequencing (DRS) data. The figure should be suitable for Nature-quality publication with clean, professional design.

---

## Layout Structure (Top to Bottom, 4 Phases)

**Overall Style:**
- Clean white or light gray background
- Modern flat design with subtle gradients
- Rounded rectangles for process boxes
- Arrows connecting phases (flowing downward)
- Professional scientific color palette
- Clear typography (sans-serif, e.g., Helvetica/Arial)

---

## Phase 1: INPUT DATA (Top Section)

**Position:** Top of figure
**Background Color:** Light blue (#E3F2FD)
**Border:** Blue (#1565C0)

**Components (arranged horizontally):**
1. **FASTQ Icon** - Document icon with "FASTQ" label, subtitle: "Raw ONT reads"
2. **BLOW5 Icon** - Document icon with "BLOW5" label, subtitle: "Raw signal data"
3. **Reference Files** - Double document icon, labels: "ref.fa", "ref.gtf"
4. **samples.tsv** - Table icon, subtitle: "Native + Control"
5. **Truth Set** - Star/bullseye icon, label: "Ground Truth", **HIGHLIGHT with golden border**

---

## Phase 2: DATA PREPARATION (Second Section)

**Position:** Below Input
**Background Color:** Light orange (#FFF3E0)
**Border:** Orange (#E65100)
**Label:** "Phase 1: Data Preparation"

**Process Flow (left to right):**
```
[FASTQ] → [minimap2] → [Filter] → [f5c index] → [f5c eventalign]
                        ↓
                   [BAM files]
```

**Key Elements:**
- Minimap2 box with alignment visualization icon
- f5c eventalign box with signal waveform icon
- **HIGHLIGHT:** f5c eventalign (critical step, use brighter color)

---

## Phase 3: MODIFICATION DETECTION (Third Section - LARGEST)

**Position:** Center of figure
**Background Color:** Light purple (#F3E5F5)
**Border:** Purple (#4A148C)
**Label:** "Phase 2: Modification Detection"

### Two Parallel Sub-sections:

**Sub-section A: Per-Comparison Tools (Left)**
- **Label:** "Comparative Tools (Native vs Control)"
- **Border:** Dashed line to indicate comparison requirement
- Tools (arranged in 2 columns, 4 rows):
  - xPore, Nanocompore
  - Baleen, pyBaleen
  - DiffErr, DRUMMER
  - ELIGOS2, EpiNano
  - psipore

**Sub-section B: Per-Sample Tools (Right)**
- **Label:** "Single-Sample Tools (ML/DL)"
- **Border:** Solid line
- **HIGHLIGHT:** Different background shade
- Tools (arranged in 2 columns):
  - TandemMod, DirectRM
  - m6ATM, Rnano
  - NanoPSU, NanoMUD
  - Penguin

**Visual Distinction:**
- Comparative tools: Blue/green tones
- Single-sample tools: Orange/red tones (ML/DL emphasis)

---

## Phase 4: POST-PROCESSING (Fourth Section)

**Position:** Below Detection
**Background Color:** Light teal (#E0F2F1)
**Border:** Teal (#00695C)
**Label:** "Phase 3: Post-Processing"

**Components:**
1. **Format Normalization** - Box with transformation icon
2. **GTF Annotation** - Box with annotation icon
3. **Output:** `{tool}_results.tsv` (document icon)

---

## Phase 5: BENCHMARKING (Bottom Section - HIGHLIGHTED)

**Position:** Bottom of figure
**Background Color:** Light green (#E8F5E9)
**Border:** Green (#1B5E20), **THICK BORDER for emphasis**
**Label:** "Phase 4: Benchmarking"
**★ THIS IS THE KEY INNOVATION - MAKE IT STAND OUT**

### Two Parallel Modules:

**Module A: Multi-Threshold Evaluation (Left)**
- **Icon:** Threshold/slider icon
- **Process:**
  - Score distribution analysis
  - Multiple threshold testing
  - Optimal threshold discovery (F1 maximization)
- **Output icons:** threshold_evaluation.tsv, optimal_thresholds.tsv

**Module B: Accuracy Metrics (Right)**
- **Icon:** Chart/metrics icon
- **Metrics displayed as small badges:**
  - Precision, Recall, F1
  - Specificity, MCC
  - AUROC, AUPRC
- **Output icons:** accuracy_summary.tsv

**Module C: Resource Aggregation (Bottom)**
- **Icon:** CPU/memory icon
- **Metrics:** Runtime, Memory, I/O
- **Output:** resource_summary.tsv

---

## FINAL OUTPUT (Bottom-most)

**Position:** Very bottom
**Background Color:** Light pink (#FCE4EC)
**Border:** Pink (#880E4F)

**Output Files (arranged horizontally):**
1. accuracy_summary.tsv
2. optimal_thresholds.tsv
3. threshold_evaluation.tsv
4. resource_summary.tsv
5. Visualization plots (small chart icons)

---

## Connecting Arrows

- **Main flow:** Gray arrows (#616161), medium thickness
- **Data to Detection:** Split arrows to both tool categories
- **Detection to Post-processing:** Converging arrows
- **Post-processing to Benchmarking:** Bold green arrow (**HIGHLIGHT**)
- **Truth Set to Benchmarking:** Golden dashed arrow (**HIGHLIGHT** - this is critical connection)

---

## Color Legend (Bottom Right Corner)

Small legend box:
- 🔵 Input Data (Blue)
- 🟠 Data Preparation (Orange)
- 🟣 Modification Detection (Purple)
- 🟢 Benchmarking (Green) ← **Key Module**
- 🔴 Output (Pink)

---

## Key Highlights for AI Generation

**CRITICAL ELEMENTS TO EMPHASIZE:**

1. **Benchmarking Module** - Should be visually prominent
   - Larger size than other phases
   - Thicker border
   - Slight glow/shadow effect

2. **Multi-Threshold Evaluation** - Unique feature
   - Add small graph showing threshold curve
   - Label: "Optimal Threshold Discovery"

3. **Truth Set Connection** - Golden dashed line
   - Connects directly to benchmarking
   - Star icon for ground truth

4. **Tool Categories** - Clear visual separation
   - Comparative vs Single-sample
   - Different color schemes

---

## Typography

- **Title:** "NanoRNAMod Benchmarking Workflow" (Bold, 24pt)
- **Phase labels:** Semi-bold, 14pt
- **Tool names:** Regular, 10pt
- **Metric names:** Italic, 9pt

---

## Additional Scientific Elements

- Add small RNA helix icon near title
- Add nanopore sequencing icon in input section
- Add DNA/RNA modification symbols (m6A, Ψ) as decorative elements
- Add small ROC curve thumbnail in benchmarking section

---

## Figure Dimensions

- Recommended: 11 inches × 14 inches (portrait)
- Resolution: 300 DPI minimum for print
- Format: Vector (SVG/PDF) preferred, or high-res PNG

---

## Prompt Examples for AI Image Generation

### For DALL-E / ChatGPT:
```
Create a scientific workflow diagram for a bioinformatics pipeline called "NanoRNAMod Benchmarking Workflow". The diagram should be suitable for Nature journal publication.

Layout: 4 phases from top to bottom
1. Input (blue): FASTQ, BLOW5, Reference files, Truth Set
2. Data Preparation (orange): minimap2 → filtering → f5c eventalign
3. Modification Detection (purple): Two parallel sections - Comparative tools (xPore, Nanocompore, Baleen, etc.) and Single-sample ML tools (TandemMod, m6ATM, etc.)
4. Benchmarking (green, HIGHLIGHTED): Multi-threshold evaluation, Accuracy metrics (Precision, Recall, F1, MCC, AUROC, AUPRC), Resource aggregation

Style: Clean, modern flat design with rounded rectangles. Professional scientific color palette. White background. Arrows connecting phases downward. The Benchmarking section should be visually emphasized as the key innovation.

Include: RNA helix decoration, nanopore icon, small ROC curve in benchmarking section, golden dashed line from Truth Set to Benchmarking.
```

### For Midjourney:
```
Scientific workflow diagram, bioinformatics pipeline, NanoRNAMod Benchmarking, Nature journal quality, clean flat design, 4 phases top to bottom, blue orange purple green color scheme, RNA modification detection, nanopore sequencing, professional scientific illustration, vector style, white background, --ar 11:14 --v 6
```

### For Stable Diffusion:
```
Professional scientific workflow diagram, bioinformatics, RNA modification detection benchmarking, Oxford Nanopore sequencing, clean modern design, 4-phase pipeline, publication quality, flat design, rounded rectangles, connecting arrows, professional color palette, high resolution, 300 DPI
```
