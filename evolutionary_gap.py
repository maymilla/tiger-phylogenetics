import numpy as np


def get_branch_length(edges, node_name):
    for parent, child, branch_length in edges:
        if child == node_name:
            return branch_length
    return None


def compute_evolutionary_loss(edges, extinct_species):
    total_tree_length = sum(branch for _, _, branch in edges)
    
    loss_per_species = {}
    for sp in extinct_species:
        bl = get_branch_length(edges, sp)
        if bl is not None:
            loss_per_species[sp] = bl
        else:
            print(f"Warning: {sp} tidak ditemukan di pohon")
            loss_per_species[sp] = 0.0
    
    total_loss = sum(loss_per_species.values())
    proportion_loss = (total_loss / total_tree_length * 100) if total_tree_length > 0 else 0.0
    
    return {
        "loss_per_species": loss_per_species,
        "total_loss": total_loss,
        "total_tree_length": total_tree_length,
        "proportion_loss": proportion_loss
    }


def compute_evolutionary_distinctiveness(matrix, species):
    scores = {}
    for i, name in enumerate(species):
        row = [matrix[i][j] for j in range(len(species)) if i != j]
        scores[name] = np.mean(row)
    
    scores = dict(sorted(scores.items(), key=lambda x: x[1], reverse=True))
    return scores


def summarize_gap_analysis(evolutionary_loss, distinctiveness_scores, node_support):
    print("\n" + "="*60)
    print("EVOLUTIONARY GAP ANALYSIS")
    print("="*60)
    
    print("\n[1] Evolutionary Loss akibat Kepunahan")
    for sp, loss in evolutionary_loss["loss_per_species"].items():
        short_name = sp.split("|")[0].strip()
        print(f"  {short_name}: {loss:.6f} substitusi/situs")

    print(f"  Total loss   : {evolutionary_loss['total_loss']:.6f}")
    print(f"  Total pohon  : {evolutionary_loss['total_tree_length']:.6f}")
    print(f"  Proporsi loss: {evolutionary_loss['proportion_loss']:.2f}%")
    
    # interpretation
    if evolutionary_loss["proportion_loss"] > 5:
        print("-> Kepunahan menyebabkan gap evolusioner SIGNIFIKAN")
    else:
        print("-> Gap evolusioner relatif kecil")
    
    print("\n[2] Evolutionary Distinctiveness (tinggi = lebih unik)")
    for sp, score in distinctiveness_scores.items():
        short_name = sp.split("|")[0].strip()
        print(f"  {short_name}: {score:.6f}")
    
    top_species = list(distinctiveness_scores.keys())[0].split("|")[0].strip()
    print(f"-> {top_species} adalah subspesies paling unik secara genetik")
    
    print("\n[3] Bootstrap Support tiap Node")
    for node, support in node_support.items():
        status = "kuat" if support >= 70 else "lemah"
        print(f"  {node}: {support}% ({status})")
    
    print("="*60)
    
    return {
        "evolutionary_loss": evolutionary_loss,
        "distinctiveness": distinctiveness_scores,
        "bootstrap_support": node_support
    }

def main(edges, matrix, species, node_support={}):
    extinct = [sp for sp in species if "sondaica" in sp.lower()]
    if not extinct:
        print("Warning: tidak ada spesies punah ditemukan di data")
        return None, None
    
    print(f"Spesies punah: {[s.split('|')[0].strip() for s in extinct]}")
    loss = compute_evolutionary_loss(edges, extinct)
    distinctiveness = compute_evolutionary_distinctiveness(matrix, species)
    summarize_gap_analysis(loss, distinctiveness, node_support)
    
    return loss, distinctiveness