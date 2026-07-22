# ACRC Report Design System — Spec

University of Illinois Air Conditioning and Refrigeration Center. White report pages, navy structure, rare orange accents, Arial throughout.

## Colors

### Brand
- Illini Blue / navy: `#13294B`
- Navy band (deck header): `#1B4284`
- Navy deep: `#0E1F3A`
- Illini Orange: `#FF5F05`
- Orange deep: `#D94F00`

### Chart series (data-viz only — never text/fills)
- Series 1 (navy): `#13294B`
- Series 2 (Storm): `#707372`
- Series 3 (Altgeld): `#C84113`

### Grayscale
- Black: `#000000`
- Ink strong: `#1F1F1F`
- Ink (body text): `#232323`
- Grey 1 (lightest fill / zebra): `#F3F3F3`
- Grey 2 (panel fills, dividers): `#EBEBEB`
- Grey 3 (hairline borders): `#CECECE`
- Grey 4 (muted labels, disabled): `#9C9C9C`
- Grey 5 (secondary text, min. grey for text): `#676767`
- White: `#FFFFFF`

### Semantic
- Text body: `#232323` · strong: `#1F1F1F` · secondary: `#676767` · muted: `#9C9C9C`
- Text on navy / on orange: `#FFFFFF`
- Link: `#13294B`
- Surface page / panel: `#FFFFFF` · soft: `#F3F3F3` · fill: `#EBEBEB` · header/dark: `#13294B`
- Border hairline: `#CECECE` · soft: `#EBEBEB` · strong: `#9C9C9C`
- Brand: `#13294B` · Accent: `#FF5F05` (sparingly) · Accent deep: `#D94F00`

## Typography

Single family: `Arial, Helvetica, 'Helvetica Neue', sans-serif`. No webfonts. Hierarchy by weight and size only.

### Web/specimen scale
- Display: 40px / 700 / lh 1.06 / ls −0.4px / `#1F1F1F`
- Title: 28px / 700 / lh 1.12 / ls −0.2px / `#1F1F1F`
- Heading: 20px / 700 / lh 1.25
- Eyebrow: 13px / 600 / lh 1.2 / ls 0.08em / uppercase
- Body-lg: 18px / 400 / lh 1.5
- Body: 16px / 400 / lh 1.55
- Caption: 13px / 400 / lh 1.4 / `#232323`
- Footer: 11px / 400 / lh 1.3 / `#676767`
- Identity line: 20px / 400 / lh 1.35; large: 30px / 400 / lh 1.2
- Numerics: `font-variant-numeric: tabular-nums`

### Slide scale (1280×720)
- Band eyebrow (section title): 38px / 700 / ls 0.2px / `#FFFFFF` / centered
- Subtitle: 26px / 700 / `#1F1F1F`
- Lede: 22px / 400 / lh 1.5 / max-width 62ch
- Bullets: 22px / lh 1.4
- Cover title: 40px / 700 / lh 1.12 / ls −0.4px / max-width 24ch
- Cover center-name: 30px / 400 / lh 1.2 / `#13294B`
- Cover uni line: 12px / 600 / ls 0.16em / uppercase / `#676767`
- Cover meta: 16px (keys 700)
- Figure caption: 18px / lh 1.35 / centered
- Sensor name: 15px / 700 · sensor sub: 12px
- Ack name: 16px / 700 · role: 13px / `#676767`
- Footer: 11px / `#676767` / centered · Page number: 14px / `#676767` / tabular-nums

## Spacing (8px base)
- xxs 2px · xs 4px · sm 8px · md 12px · lg 16px · xl 24px · xxl 32px · huge 64px
- Container max: 1100px
- Slide content inset: 32px
- Header band height (web): 72px
- Section padding-y (web): 64px
- Card padding: 24px

## Border radius
- xs 2px · sm 4px (chips) · md 6px (cards, figure frames) · lg 8px (panels) · xl 10px · pill 9999px (tags)

## Elevation
- Shadow 1: `0 1px 2px rgba(0,0,0,0.08)`
- Shadow 2: `0 4px 16px rgba(0,0,0,0.10)`
- Default: flat; hairline borders over shadows.

## Slide layout rules (1280×720)

### Page
- White background, `overflow: hidden`, Arial, ink `#232323`.

### Navy header band
- Absolute top, full width, height 90px, background `#13294B`, padding 0 40px.
- Centered white 38px/700 section title.
- Block-I logo absolute right 40px, height 61px, vertically centered.

### Content stage
- Absolute below band: top 90px to bottom, padding 26px 40px 54px, flex column.

### Bullets
- No list markers; 20px gap between items; 28px left padding.
- Square navy bullet: 11×11px `#13294B`, at left 0 / top 10px.

### Accomplishments grid (bullets + figure)
- 2-col grid `42% / 1fr`, gap 48px, full height; bullets gap 26px, self-centered; figure `object-fit: contain`, 18px caption 12px below, figure number bold navy.

### Figure rows
- Flex row, gap 24px, margin-top 16px; each figure flex column, frame fills with `object-fit: contain`; caption 18px centered, 10px above margin.

### Sensor strip
- 4-col grid, gap 18px; frames square (aspect-ratio 1/1, `object-fit: cover`); name 15px/700 8px below, sub 12px, centered.

### Cover slide
- Centered flex column, padding 0 90px.
- Block-I 66px → uni line → center-name, head group gap 8px.
- Orange rule: 120×3px `#FF5F05`, margin 22px 0.
- Title 40px/700, then meta stack margin-top 26px, gap 6px.
- Headshots: top corners at 40px inset, 150px wide columns; avatar 96px circle, 2px `#13294B` border, `#EBEBEB` fill; name 13px.

### Acknowledgements
- 3-col grid, gap 22px 28px, margin-top 20px.
- Avatar 84px circle, 2px `#CECECE` border, `#EBEBEB` fill; name 16px/700, role 13px `#676767`.

### Footer / page number
- Footer: absolute bottom 16px, centered, 11px `#676767`.
- Page number: absolute right 40px / bottom 14px, 14px `#676767`, tabular-nums.

### Charts
- Axis stroke `#CECECE` 1px; gridline stroke `#EBEBEB` 1px; series colors per chart scale above.
