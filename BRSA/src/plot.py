import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import mplcursors
import plotly.express as px
import plotly.graph_objects as go
import gseapy as gp
import os

try:
    from scipy.cluster.hierarchy import dendrogram
except ImportError:
    pass

try:
    from lifelines import KaplanMeierFitter
except ImportError:
    pass

try:
    import networkx as nx
except ImportError:
    pass

# Slaps some paint on the canvas to show how conditions spread out in the multidimensional void
def paint_dimensional_void_masterpiece(ax, dim_data, metadata_df):
    conditions = metadata_df['condition'].values
    unique_conditions = np.unique(conditions)
    colors = plt.cm.tab10(np.linspace(0, 1, len(unique_conditions)))
    color_map = dict(zip(unique_conditions, colors))
    
    for i, condition in enumerate(conditions):
        ax.scatter(
            dim_data['components'][i, 0],
            dim_data['components'][i, 1],
            color=color_map[condition],
            label=condition if condition not in ax.get_legend_handles_labels()[1] else "",
            alpha=0.8,
            s=80, edgecolors='w'
        )
        
    dtype = dim_data['type']
    if dtype == "PCA":
        ax.set_xlabel(f"PC1 ({dim_data['explained_variance'][0]*100:.1f}%)", fontweight='bold')
        ax.set_ylabel(f"PC2 ({dim_data['explained_variance'][1]*100:.1f}%)", fontweight='bold')
    else:
        ax.set_xlabel(f"{dtype} Dimension 1", fontweight='bold')
        ax.set_ylabel(f"{dtype} Dimension 2", fontweight='bold')
        
    ax.set_title(f"{dtype} Projection", fontweight='bold', fontsize=14)
    ax.legend(title="Conditions")
    ax.grid(True, linestyle='--', alpha=0.6)

# Summons the true hierarchical demon with dendrograms for samples and genes
def summon_clustermap_demon(heatmap_matrix, metadata_df, color_palette="mako"):
    if 'sample' in metadata_df.columns:
        condition_series = metadata_df.set_index('sample')['condition']
    else:
        condition_series = metadata_df['condition']
        
    lut = dict(zip(condition_series.unique(), sns.color_palette("Set2", len(condition_series.unique()))))
    col_colors = heatmap_matrix.columns.map(condition_series).map(lut)
    
    g = sns.clustermap(
        heatmap_matrix,
        cmap=color_palette,
        col_colors=col_colors,
        z_score=0, 
        figsize=(10, 8),
        xticklabels=True,
        yticklabels=False,
        linewidths=0
    )
    g.fig.suptitle("Hierarchical Clustermap (Top Variable Genes)", fontweight='bold', fontsize=14, y=1.05)
    return g.fig

# Blows up the significant genes to the top like a glorious volcano eruption
def erupt_volcano_plot(ax, de_results, pval_cutoff):
    de_results['neg_log10_pval'] = -np.log10(de_results['p-adj'] + 1e-300)
    colors = ['#444444' if not sig else '#e74c3c' for sig in de_results['significant']]
    
    sc = ax.scatter(
        de_results['log2fc'],
        de_results['neg_log10_pval'],
        c=colors,
        alpha=0.6,
        s=25,
        edgecolors='none'
    )
    
    ax.axhline(-np.log10(pval_cutoff), linestyle='--', color='#2c3e50', alpha=0.8)
    ax.axvline(0, linestyle='--', color='#2c3e50', alpha=0.8)
    ax.set_xlabel("Log2 Fold Change", fontweight='bold')
    ax.set_ylabel("-log10 FDR (p-adj)", fontweight='bold')
    ax.set_title("Volcano Plot", fontweight='bold', fontsize=14)
    
    cursor = mplcursors.cursor(sc, hover=True)
    @cursor.connect("add")
    def on_add(sel):
        gene = de_results.index[sel.target.index]
        sel.annotation.set_text(
            f"Gene: {gene}\n"
            f"LFC: {de_results.loc[gene, 'log2fc']:.2f}\n"
            f"FDR: {de_results.loc[gene, 'p-adj']:.2e}"
        )
        sel.annotation.get_bbox_patch().set_alpha(0.8)

# Paints the constellation of expression vs fold change to judge the normalization gods
def paint_ma_constellation(ax, de_results):
    colors = ['#444444' if not sig else '#3498db' for sig in de_results['significant']]
    
    base_mean = de_results.get('mean_exp', (de_results['mean_1'] + de_results['mean_2']) / 2)
    log_base_mean = np.log2(base_mean + 1)
    
    sc = ax.scatter(
        log_base_mean,
        de_results['log2fc'],
        c=colors,
        alpha=0.6,
        s=20,
        edgecolors='none'
    )
    
    ax.axhline(0, linestyle='--', color='#e74c3c', alpha=0.8)
    ax.set_xlabel("Log2 Mean Expression", fontweight='bold')
    ax.set_ylabel("Log2 Fold Change", fontweight='bold')
    ax.set_title("MA Constellation Plot", fontweight='bold', fontsize=14)
    
    cursor = mplcursors.cursor(sc, hover=True)
    @cursor.connect("add")
    def on_add(sel):
        gene = de_results.index[sel.target.index]
        sel.annotation.set_text(f"Gene: {gene}\nLFC: {de_results.loc[gene, 'log2fc']:.2f}")

# Visualizes the WGCNA topological overlap dendrogram
def plot_dendrogram_of_madness(ax, linkage_matrix):
    dendrogram(linkage_matrix, ax=ax, no_labels=True, color_threshold=0)
    ax.set_title("WGCNA Co-expression Dark Tapestry", fontweight='bold', fontsize=14)
    ax.set_ylabel("Distance")
    ax.set_xlabel("Co-expressed Souls (Genes)")

# Draws the grim curve of survival based on a specific gene's expression
def predict_mortal_doom_curves(ax, metadata, expr_matrix, gene_name, time_col='time', event_col='event'):
    try:
        kmf = KaplanMeierFitter()
        gene_expr = expr_matrix.loc[gene_name]
        median_expr = gene_expr.median()
        
        high_expr_samples = gene_expr[gene_expr > median_expr].index
        low_expr_samples = gene_expr[gene_expr <= median_expr].index
        
        T = metadata[time_col]
        E = metadata[event_col]
        
        kmf.fit(T[high_expr_samples], event_observed=E[high_expr_samples], label=f"High {gene_name} Curse")
        kmf.plot_survival_function(ax=ax, color='#e74c3c')
        
        kmf.fit(T[low_expr_samples], event_observed=E[low_expr_samples], label=f"Low {gene_name} Blessing")
        kmf.plot_survival_function(ax=ax, color='#3498db')
        
        ax.set_title(f"Mortal Doom Curve (Kaplan-Meier): {gene_name}", fontweight='bold', fontsize=14)
        ax.set_ylabel("Survival Probability")
        ax.set_xlabel("Time Runic Unit")
    except Exception as e:
        ax.text(0.5, 0.5, f"Cannot predict Doom Curve:\n{e}", ha='center', va='center', color='white')

# Interrogates the spirits of Random Forest feature importances
def interrogate_feature_importance_spirits(ax, importance_series):
    importance_series.sort_values(ascending=True).plot(kind='barh', ax=ax, color='#8a2be2', edgecolor='#292e42')
    ax.set_title("ML Golems: Top Predictive Genes", fontweight='bold', fontsize=14)
    ax.set_xlabel("Feature Importance Score")

# Visualizes the chaotic hairball of protein interactions
def paint_hairball_of_chaos(ax, ppi_graph):
    if len(ppi_graph.nodes) == 0:
        ax.text(0.5, 0.5, "The Interactome Spirits are silent (No network data).", ha='center', va='center', color='white')
        return
        
    pos = nx.spring_layout(ppi_graph, seed=42)
    nx.draw_networkx_nodes(ppi_graph, pos, ax=ax, node_color='#ff8c00', node_size=100, alpha=0.8)
    nx.draw_networkx_edges(ppi_graph, pos, ax=ax, alpha=0.5, edge_color='#a9b1d6')
    nx.draw_networkx_labels(ppi_graph, pos, ax=ax, font_size=8, font_color='white')
    ax.set_title("STRING DB: Hairball of Chaos", fontweight='bold', fontsize=14)
    ax.axis('off')

# Uncovers hidden isoform variances
def uncover_hidden_isoform_mutations(ax, frankenstein_isoforms, expr_matrix):
    data = expr_matrix.loc[frankenstein_isoforms[:20]].T
    sns.boxplot(data=data, ax=ax, orient='h', palette='magma')
    ax.set_title("Frankenstein Limbs: Isoform Expression Variance", fontweight='bold', fontsize=14)
    ax.set_xlabel("Normalized Expression")

# Renders the beautiful running enrichment score mountain
def paint_gsea_mountain(ax, gsea_obj, term_name):
    try:
        gp.plot.gseaplot(rank_metric=gsea_obj.ranking, term=term_name, **gsea_obj.results[term_name], ax=ax)
    except Exception as e:
        ax.text(0.5, 0.5, f"Cannot render GSEA mountain:\n{e}", ha='center', va='center', color='white')

# Forges standalone interactive HTML files using Plotly
def forge_interactive_plotly_artifacts(treasures: dict, output_dir: str):
    os.makedirs(output_dir, exist_ok=True)
    
    # 1. Interactive Dimensionality Reduction
    dim_data = treasures['dim_red']
    df_dim = pd.DataFrame(dim_data['components'], columns=['Dim1', 'Dim2'])
    
    if 'sample' in treasures['metadata'].columns:
        df_dim['Sample'] = treasures['metadata']['sample'].values
        df_dim['Condition'] = treasures['metadata']['condition'].values
    else:
        df_dim['Sample'] = treasures['metadata'].index.values
        df_dim['Condition'] = treasures['metadata']['condition'].values
    
    fig_dim = px.scatter(df_dim, x='Dim1', y='Dim2', color='Condition', hover_name='Sample', 
                         title=f"Interactive {dim_data['type']} Projection", template="plotly_dark")
    fig_dim.write_html(os.path.join(output_dir, f"interactive_{dim_data['type'].lower()}.html"))
    
    # 2. Interactive Volcano
    de = treasures['de_results'].reset_index()
    de['neg_log10_pval'] = -np.log10(de['p-adj'] + 1e-300)
    de['Significance'] = ['Significant' if s else 'Not Significant' for s in de['significant']]
    
    fig_volc = px.scatter(de, x='log2fc', y='neg_log10_pval', color='Significance', hover_name='gene',
                          color_discrete_map={'Significant': '#e74c3c', 'Not Significant': '#7f8c8d'},
                          title="Interactive Volcano Plot", template="plotly_dark")
    fig_volc.write_html(os.path.join(output_dir, "interactive_volcano.html"))
    
    # 3. Interactive MA Constellation
    de['base_mean_log2'] = np.log2(de.get('mean_exp', (de['mean_1'] + de['mean_2']) / 2) + 1)
    fig_ma = px.scatter(de, x='base_mean_log2', y='log2fc', color='Significance', hover_name='gene',
                        color_discrete_map={'Significant': '#3498db', 'Not Significant': '#7f8c8d'},
                        title="Interactive MA Constellation Plot", template="plotly_dark")
    fig_ma.write_html(os.path.join(output_dir, "interactive_ma_constellation.html"))
