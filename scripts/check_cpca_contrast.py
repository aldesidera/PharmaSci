from __future__ import annotations


def channel(value: int) -> float:
    value /= 255
    return value / 12.92 if value <= 0.04045 else ((value + 0.055) / 1.055) ** 2.4


def luminance(hex_color: str) -> float:
    color = hex_color.lstrip('#')
    rgb = [int(color[i:i + 2], 16) for i in (0, 2, 4)]
    return 0.2126 * channel(rgb[0]) + 0.7152 * channel(rgb[1]) + 0.0722 * channel(rgb[2])


def contrast(foreground: str, background: str) -> float:
    a, b = luminance(foreground), luminance(background)
    return (max(a, b) + 0.05) / (min(a, b) + 0.05)

pairs = [
    ('#9a3412', '#fff1eb', 'candidato laranja escuro 1'),
    ('#7c2d12', '#fff1eb', 'candidato laranja escuro 2'),
    ('#ea580c', '#fff1eb', 'texto ativo laranja sobre fundo ativo'),
    ('#3d2020', '#fff1eb', 'texto principal sobre fundo ativo'),
    ('#6f514b', '#fff1eb', 'texto secundário sobre fundo ativo'),
    ('#b45309', '#fff7ed', 'texto de revisão sobre fundo de revisão'),
    ('#3d2020', '#ffffff', 'texto principal sobre fundo branco'),
]
for foreground, background, label in pairs:
    ratio = contrast(foreground, background)
    print(f'{label}: {foreground} / {background} = {ratio:.2f}:1')
