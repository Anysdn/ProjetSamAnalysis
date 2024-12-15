import sys
from collections import defaultdict, Counter
import matplotlib.pyplot as plt


# Function to parse a SAM file and yield each read as a dictionary
def parse_sam(fichier_sam): 
    with open(fichier_sam, 'r') as file: 
        for ligne in file: 
            # Skip header lines (starting with '@')
            if ligne[0] != '@':
                colonnes = ligne.strip().split('\t')  # Split the line by tab characters
                yield {
                    "QNAME": colonnes[0],   # Read name
                    "FLAG": int(colonnes[1]),  # Alignment status flag
                    "RNAME": colonnes[2],  # Reference sequence name (chromosome)
                    "POS": int(colonnes[3]),  # Alignment position on the sequence
                    "MAPQ": int(colonnes[4]),  # Mapping quality
                    "CIGAR": colonnes[5],  # CIGAR string
                }

# Function to check if a read is mapped (FLAG bit 3 not set)
def verif_mappe(flag):
    return (flag & 4) == 0  # If bit 3 is 0, the read is mapped

# Function to check if a read is partially mapped (based on MAPQ and CIGAR string)
def est_partiellement_mappé(mapq, cigar):
    """
    A read is considered partially mapped if it has:
    - Low MAPQ (<30)
    - 'S' (soft clipping) or 'I'/'D' in the CIGAR string (indels).
    """
    return mapq < 30 or 'S' in cigar or 'I' in cigar or 'D' in cigar

# Function to check if the read belongs to a properly paired alignment (FLAG bits 1 and 2 set)
def est_paire_correcte(flag):
    return (flag & 1) != 0 and (flag & 2) != 0  # Both bits must be set for proper pairing

# Main function to count reads and compute statistics
def comptage_reads(fichier_sam):
    stats = {
        "total_reads": 0,  # Total number of reads
        "mapped_reads": 0,  # Reads that are mapped
        "unmapped_reads": 0,  # Reads that are unmapped
        "paired_reads": 0,  # Reads that are paired
        "partially_mapped_reads": 0,  # Reads partially mapped
        "reads_par_chromosome": Counter(),  # Count reads by chromosome
        "mapq_distribution": Counter(),  # Distribution of MAPQ scores
        "cigar_distribution": Counter(),  # Distribution of CIGAR alignment types
    }
    
    paired_groups = defaultdict(list)  # Group reads by their QNAME (read name)

    for read in parse_sam(fichier_sam): 
        stats["total_reads"] += 1  # Increment total read count
        flag = read["FLAG"]
        rname = read["RNAME"]
        mapq = read["MAPQ"]
        cigar = read["CIGAR"]

        # Analyze the CIGAR string for alignment type
        type_alignement = "non_match"
        match_found = False  # Indicates if a "simple match" is found
        for i in range(len(cigar)):
            if cigar[i].isdigit():
                j = i
                while j < len(cigar) and cigar[j].isdigit():
                    j += 1
                if j < len(cigar) and cigar[j] == 'M':
                    match_found = True
        if match_found:
            type_alignement = "simple_match"

        stats["cigar_distribution"][type_alignement] += 1

        # Check if the read is mapped
        if verif_mappe(flag): 
            stats["mapped_reads"] += 1
            stats["reads_par_chromosome"][rname] += 1
            if est_partiellement_mappé(mapq, cigar): 
                stats["partially_mapped_reads"] += 1
        else:
            stats["unmapped_reads"] += 1 

        # Group reads by QNAME
        paired_groups[read["QNAME"]].append(read)

        # Group MAPQ scores into bins of 10
        tranche = (mapq // 10) * 10
        stats["mapq_distribution"][f"[{tranche}-{tranche + 9}]"] += 1

    return stats, paired_groups 

# Function to analyze paired reads
def analyse_paires(paired_groups):
    pair_stats = {
        "paired_reads": 0,  # Total number of read pairs
        "pairs_one_mapped_one_unmapped": 0,  # Pairs with one mapped and one unmapped read
        "pairs_one_mapped_one_partially_mapped": 0,  # Pairs with one mapped and one partially mapped read
    }    
    
    for qname, group in paired_groups.items():
        if len(group) != 2:  # Skip groups that don't contain exactly two reads
           continue

        pair_stats["paired_reads"] += 1
        read1, read2 = group
        flag1, flag2 = read1["FLAG"], read2["FLAG"]
        mapq1, mapq2 = read1["MAPQ"], read2["MAPQ"]
        cigar1, cigar2 = read1["CIGAR"], read2["CIGAR"]

        # Check for various pairing conditions
        if verif_mappe(flag1) and not verif_mappe(flag2):
            pair_stats["pairs_one_mapped_one_unmapped"] += 1
        elif not verif_mappe(flag1) and verif_mappe(flag2):
            pair_stats["pairs_one_mapped_one_unmapped"] += 1
        elif verif_mappe(flag1) and est_partiellement_mappé(mapq2, cigar2):
            pair_stats["pairs_one_mapped_one_partially_mapped"] += 1
        elif est_partiellement_mappé(mapq1, cigar1) and verif_mappe(flag2):
            pair_stats["pairs_one_mapped_one_partially_mapped"] += 1

    return pair_stats 

# Function to display the results
def afficher_resultats(stats, pair_stats):
    print("\n=== Results of Analysis ===")
    print(f"Total reads: {stats['total_reads']}")
    print(f"Mapped reads: {stats['mapped_reads']}")
    print(f"Unmapped reads: {stats['unmapped_reads']}")
    print(f"Paired reads: {pair_stats['paired_reads']}")
    print(f"Partially mapped reads: {stats['partially_mapped_reads']}")
    
    print("\n=== Pair Analysis ===")
    print(f"Pairs with one mapped and one unmapped: {pair_stats['pairs_one_mapped_one_unmapped']}")
    print(f"Pairs with one mapped and one partially mapped: {pair_stats['pairs_one_mapped_one_partially_mapped']}")

    print("\nReads by chromosome:")
    for chromosome, count in stats["reads_par_chromosome"].items():
        print(f"  {chromosome}: {count} reads")
    
    print("\nMAPQ score distribution:")
    for mapq, count in sorted(stats["mapq_distribution"].items()):
        print(f"  MAPQ {mapq}: {count} reads")
    
    print("\nCIGAR alignment types:")
    for alignement, count in stats["cigar_distribution"].items():
        print(f"  {alignement}: {count} reads")


def tracer_distribution_mapq(mapq_distribution):
    # Trier les bins de MAPQ pour l'affichage correct
    bins = sorted(mapq_distribution.keys())
    valeurs = [mapq_distribution[bin] for bin in bins]
    
    # Création du graphique
    plt.figure(figsize=(10, 6))
    plt.bar(bins, valeurs, color='skyblue', edgecolor='black')
    plt.xlabel('Tranches de MAPQ', fontsize=12)
    plt.ylabel('Nombre de lectures', fontsize=12)
    plt.title('Distribution des scores MAPQ', fontsize=14)
    plt.xticks(rotation=45, fontsize=10)
    plt.tight_layout()
    
    # Sauvegarder ou afficher le graphique
    plt.savefig("mapq_distribution.png")  # Sauvegarde du graphique sous forme d'image
    plt.show()  # Affiche le graphique à l'écran


if __name__ == "__main__":
    # Check command-line arguments
    if len(sys.argv) < 2:
        print("Specify a SAM file as an argument!")
        sys.exit(1)

    fichier_sam = sys.argv[1]
    
    stats, paired_groups = comptage_reads(fichier_sam)
    pair_stats = analyse_paires(paired_groups)
    afficher_resultats(stats, pair_stats)

    tracer_distribution_mapq(stats["mapq_distribution"])
