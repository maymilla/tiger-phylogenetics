import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from collections import defaultdict
import os

def edges_to_tree(edges):
    adjacency = defaultdict(list)
    all_children = set()
    all_parents = set()

    for parent, child, branch_length in edges:
        adjacency[parent].append((child, branch_length))
        all_children.add(child)
        all_parents.add(parent)

    roots = all_parents - all_children
    root = list(roots)[0] if roots else list(all_parents)[0]
    leaf_nodes = all_children - all_parents

    return adjacency, root, leaf_nodes

def assign_positions(root, adjacency, leaf_nodes):
    x_positions = {}
    y_positions = {}
    leaf_counter = [0]

    def get_y(node):
        if node in leaf_nodes:
            y = leaf_counter[0]
            leaf_counter[0] += 1
            y_positions[node] = y
            return y
        child_ys = [get_y(child) for child, _ in adjacency[node]]
        y = sum(child_ys) / len(child_ys)
        y_positions[node] = y
        return y

    def get_x(node, x_current):
        x_positions[node] = x_current
        for child, branch_length in adjacency[node]:
            get_x(child, x_current + branch_length)

    get_y(root)
    get_x(root, 0)

    return {node: (x_positions[node], y_positions[node]) for node in x_positions}

def get_species_type(name):
    name_lower = name.lower()
    if "sumatrae" in name_lower:
        return "sumatrae"
    elif "sondaica" in name_lower:
        return "extinct"
    else:
        return "other"

def shorten_name(name):
    if "|" not in name:
        return name
    parts = name.split("|")
    species = parts[0].strip().replace("_", " ")
    accession = parts[2].strip() if len(parts) >= 3 else ""
    return f"{species} ({accession})" if accession else species

def draw_tree(edges, node_support=None, title="Pohon Filogenetik Harimau Nusantara dan Asia"):
    adjacency, root, leaf_nodes = edges_to_tree(edges)
    positions = assign_positions(root, adjacency, leaf_nodes)

    color_map = {
        "sumatrae": "#1D9E75",
        "extinct":  "#E24B4A",
        "other":    "#888780",
    }

    fig, ax = plt.subplots(figsize=(16, 9))

    for parent, child, _ in edges:
        if parent not in positions or child not in positions:
            continue
        x_parent, y_parent = positions[parent]
        x_child, y_child = positions[child]
        ax.plot([x_parent, x_child], [y_child, y_child],
                color="#444444", linewidth=1.2, zorder=1)
        ax.plot([x_parent, x_parent], [y_parent, y_child],
                color="#444444", linewidth=1.2, zorder=1)

    for node, (x, y) in positions.items():
        if node in leaf_nodes:
            sp_type = get_species_type(node)
            color = color_map[sp_type]
            ax.scatter(x, y, color=color, s=60, zorder=3)
            label = shorten_name(node)
            ax.text(x + 0.0003, y, label,
                    va="center", ha="left", fontsize=9,
                    color=color,
                    fontweight="bold" if sp_type in ("sumatrae", "extinct") else "normal")
        else:
            ax.scatter(x, y, color="#B4B2A9", s=25, zorder=3)
            if node_support and node in node_support:
                support = node_support[node]
                color_support = "#185FA5" if support >= 70 else "#BA7517"
                ax.text(x, y + 0.12, f"{support:.0f}%",
                        va="bottom", ha="center",
                        fontsize=7, color=color_support)

    x_max = max(x for x, y in positions.values())
    ax.set_xlim(left=-x_max * 0.05, right=x_max * 1.5)

    legend_elements = [
        mpatches.Patch(color="#1D9E75", label="Harimau Sumatera"),
        mpatches.Patch(color="#E24B4A", label="Harimau Jawa (punah)"),
        mpatches.Patch(color="#888780", label="Subspesies Asia lain"),
    ]
    ax.legend(handles=legend_elements,
              loc="upper right",
              fontsize=7,
              framealpha=0.8,
              edgecolor="#D3D1C7",
              borderpad=0.5,
              handlelength=1.0,
              handleheight=0.8)

    ax.set_title(title, fontsize=13, fontweight="bold", pad=15)
    ax.text(0.5, 1.01,
            "Angka di percabangan = bootstrap support (persentase dari 1000 ulangan); "
            "biru ≥70% (kuat), kuning <70% (lemah)",
            transform=ax.transAxes,
            ha="center", fontsize=7.5, color="#5F5E5A")
    ax.set_xlabel("Jarak Genetik K2P (substitusi per situs)", fontsize=10)
    ax.set_yticks([])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)

    plt.tight_layout()

    os.makedirs("results", exist_ok=True)
    plt.savefig("results/phylogenetic_tree.png", dpi=150, bbox_inches="tight")
    plt.show()
    print("Pohon disimpan ke results/phylogenetic_tree.png")

def main(edges, node_support=None):
    draw_tree(edges, node_support=node_support)