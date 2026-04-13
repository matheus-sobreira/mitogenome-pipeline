#!/usr/bin/env python3
"""Gera gráfico radar comparativo das ferramentas de montagem de mitogenomas."""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# ---------- eixos ----------
categories = [
    'Reprodu-\ntibilidade',
    'Completude\ndo fluxo',
    'Generalidade\ntaxonômica',
    'Automação',
    'Acessibilidade\ncomputacional',
]
N = len(categories)

# ---------- scores (escala 1-5) ----------
# Critérios definidos na Tabela tab:scores_radar do TCC
tools = {
    # Reprodutibilidade | Completude | Generalidade | Automação | Acessibilidade
    'MitoHiFi':       [2, 3, 3, 4, 1],  # Conda, montagem+detecção heteroplasmia, metazoários, bons defaults, requer HPC
    'MITObim':        [1, 1, 3, 2, 2],  # instalação manual, só montagem, metazoários, scripts isolados, servidor recomendado
    'GetOrganelle':   [2, 1, 4, 3, 2],  # Conda, só montagem, metazoários+plantas, pipeline com params, servidor recomendado
    'MToolBox':       [1, 3, 1, 3, 3],  # instalação manual, montagem+anotação+variantes, só humano, pipeline com params, desktop
    'MitoZ':          [4, 4, 3, 4, 3],  # Docker/Singularity, montagem+anotação+visualização, animais, bons defaults, desktop
    'MitoFinder':     [3, 4, 4, 3, 3],  # Singularity, montagem+anotação+GenBank, múltiplos códigos, params manuais, desktop
    'NOVOPlasty':     [1, 1, 4, 2, 3],  # instalação manual, só montagem, amplo (mito+plasto), config manual, desktop
    'Este trabalho':  [5, 5, 3, 5, 5],  # Docker+Nextflow+versões pinadas, fluxo completo, metazoários, um comando+auto, notebook 16GB
}

# Ângulos do radar
angles = np.linspace(0, 2 * np.pi, N, endpoint=False).tolist()
angles += angles[:1]  # fechar o polígono

# ---------- cores e estilos ----------
tool_styles = {
    'MitoHiFi':      {'color': '#7f8c8d', 'alpha': 0.15, 'lw': 1.5, 'ls': '-',  'marker': 'o', 'ms': 5},
    'MITObim':        {'color': '#95a5a6', 'alpha': 0.10, 'lw': 1.2, 'ls': '--', 'marker': 'o', 'ms': 4},
    'GetOrganelle':   {'color': '#bdc3c7', 'alpha': 0.10, 'lw': 1.5, 'ls': '-',  'marker': 'o', 'ms': 5},
    'MToolBox':       {'color': '#7f8c8d', 'alpha': 0.10, 'lw': 1.2, 'ls': '--', 'marker': 'o', 'ms': 4},
    'MitoZ':          {'color': '#636e72', 'alpha': 0.15, 'lw': 1.8, 'ls': '-',  'marker': 'o', 'ms': 5},
    'MitoFinder':     {'color': '#2d3436', 'alpha': 0.10, 'lw': 1.5, 'ls': '-',  'marker': 'o', 'ms': 5},
    'NOVOPlasty':     {'color': '#b2bec3', 'alpha': 0.10, 'lw': 1.2, 'ls': '--', 'marker': 'o', 'ms': 4},
    'Este trabalho':  {'color': '#e74c3c', 'alpha': 0.20, 'lw': 2.5, 'ls': '-',  'marker': 'o', 'ms': 7},
}

# ---------- plot ----------
fig, ax = plt.subplots(figsize=(10, 10), subplot_kw=dict(polar=True))
fig.patch.set_facecolor('white')

# Grade
ax.set_ylim(0, 5.5)
ax.set_yticks([1, 2, 3, 4, 5])
ax.set_yticklabels(['1', '2', '3', '4', '5'], fontsize=9, color='#888888')
ax.set_xticks(angles[:-1])
ax.set_xticklabels(categories, fontsize=12, fontweight='bold', color='#2C3E50')

# Linhas de grade
ax.yaxis.grid(True, color='#cccccc', linewidth=0.5)
ax.xaxis.grid(True, color='#cccccc', linewidth=0.5)
ax.spines['polar'].set_color('#cccccc')

# Plotar cada ferramenta
for name, scores in tools.items():
    values = scores + scores[:1]  # fechar o polígono
    style = tool_styles[name]
    ax.plot(angles, values,
            color=style['color'], linewidth=style['lw'],
            linestyle=style['ls'], marker=style['marker'],
            markersize=style['ms'], label=name, zorder=3)
    ax.fill(angles, values,
            color=style['color'], alpha=style['alpha'], zorder=2)

# Legenda
legend = ax.legend(loc='upper right', bbox_to_anchor=(1.32, 1.12),
                   fontsize=11, frameon=True, framealpha=0.9,
                   edgecolor='#cccccc')

# Título
fig.text(0.5, 0.97,
         'Comparação das ferramentas/pipelines\nem 5 eixos',
         ha='center', va='center', fontsize=16, fontweight='bold', color='#2C3E50')

plt.tight_layout(rect=[0, 0, 0.88, 0.94])
plt.savefig('/home/matheus/mitogenome-pipeline/latex/imagens/radar_comparativo.png',
            dpi=200, bbox_inches='tight', facecolor='white', pad_inches=0.3)
print("OK")
