#!/usr/bin/env python3
"""Compose Figure 2: new panel A (ML tree) with panels B and C retained from the original figure."""
import sys
from PIL import Image, ImageDraw, ImageFont
A = Image.open('figures/Figure2A.png').convert('RGB')
old = Image.open(sys.argv[1]).convert('RGB')          # original Figure 2 image
if old.width > 3000:   # high-resolution original (4088 px wide)
    ImageDraw.Draw(old).rectangle((2085, 820, 2240, 912), fill='white')   # panel-A trap icon above B
    B = old.crop((2085, 885, 3035, 1825)); C = old.crop((3060, 105, 4075, 1845))
else:                  # low-resolution copy embedded in the manuscript
    ImageDraw.Draw(old).rectangle((700, 270, 790, 318), fill='white')
    B = old.crop((705, 286, 1015, 655)); C = old.crop((1022, 22, 1377, 655))
import numpy as np
def whiten(im):
    a = np.asarray(im).copy(); m = a.min(axis=2) > 238; a[m] = 255; return Image.fromarray(a)
B, C = whiten(B), whiten(C)
H = A.height
w = int(H / (B.height / B.width + C.height / C.width))
Bs = B.resize((w, int(w * B.height / B.width)), Image.LANCZOS); Cs = C.resize((w, H - Bs.height), Image.LANCZOS)
pad = 60
canvas = Image.new('RGB', (A.width + w + pad, H), 'white')
canvas.paste(A, (0, 0)); canvas.paste(Bs, (A.width + pad, 0)); canvas.paste(Cs, (A.width + pad, Bs.height))
d = ImageDraw.Draw(canvas)
f = ImageFont.truetype('/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf', 48)
d.text((10, 10), 'A', fill='black', font=f)
d.text((A.width + 5, 5), 'B', fill='black', font=f)
d.text((A.width + 5, Bs.height + 10), 'C', fill='black', font=f)
canvas.save('figures/Figure2.png', dpi=(300, 300)); print(canvas.size)
