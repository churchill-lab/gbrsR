def plot(
        genoprob_file: str,
        output_file: str = None,
        output_format: str = 'pdf',
        sample_name: str = '',
        grid_size: int = 2,
        xt_max: int = 5000,
        xt_size: int = 500,
        grid_width: float = 0.01,
) -> None:
    """
    Generate genome reconstruction visualization plot.

    This function creates a comprehensive visualization of the GBRS genome reconstruction results,
    showing the most likely genotype at each position across all chromosomes. The plot displays
    both haplotypes for each position, allowing identification of recombination events and genomic
    structure patterns.

    The visualization uses a stacked bar chart format where each chromosome is represented by two
    horizontal bars (one for each haplotype). Different founder strains are color-coded, making it
    easy to identify regions of shared ancestry and recombination breakpoints. The plot includes
    recombination counts for each chromosome and total recombination count.

    ALGORITHM:
    - Load genotype probabilities and determine most likely genotypes
    - For each chromosome:
       - Extract most likely diplotype at each position
       - Separate into two haplotypes
       - Identify recombination events (genotype changes)
       - Create color-coded bars for visualization
    - Generate plot with proper spacing, labels, and annotations
    - Save in specified format with high resolution

    WHAT THIS DOES:
    - Visualizes GBRS genome reconstruction results
    - Shows haplotype structure across all chromosomes
    - Identifies recombination events and breakpoints
    - Provides quantitative recombination statistics
    - Enables visual comparison of genomic structure

    USE CASES:
    - GBRS visualization: Create publication-quality genome plots
    - Quality assessment: Visualize reconstruction quality and patterns
    - Recombination analysis: Identify and count recombination events
    - Population studies: Compare genomic structure across individuals
    - Publication figures: Generate high-resolution plots for papers

    Args:
        genoprob_file: Path to the genotype probability file (npz format).
            Type: str
            Description: File containing genotype probabilities from GBRS
            reconstruction or interpolation. Can be either gene-based or
            grid-based probabilities.
            Format: Compressed numpy file with one array per chromosome
            Shape: (num_diplotypes, num_positions) for each chromosome
            Example: 'DO336.genoprobs.npz' or 'DO336.interpolated.genoprobs.npz'

        output_file: Output filename for the plot.
            Type: str
            Default: None (auto-generated from input filename)
            Description: Name of the output file for the generated plot.
            If None, generates name based on input filename with 'plotted' prefix.
            Example: 'DO336.plotted.genome.pdf'

        output_format: File format for the output plot.
            Type: str
            Default: 'pdf'
            Description: Image format for the output file. Common formats
            include 'pdf', 'png', 'svg', 'jpg', 'tiff'.
            Example: 'pdf' for publication-quality vector graphics

        sample_name: Name of the sample for plot title.
            Type: str
            Default: ''
            Description: Sample identifier to include in the plot title.
            If empty, uses the filename as the sample name.
            Example: 'DO336' or 'Mouse_001'

        grid_size: Size parameter for grid scaling.
            Type: int
            Default: 2
            Description: Advanced parameter affecting the scaling of the
            x-axis grid. Used to adjust the visual spacing of positions.
            Range: Typically 1 to 5
            Example: 2 for standard spacing

        xt_max: Maximum x-axis tick value.
            Type: int
            Default: 5000
            Description: Maximum value for x-axis ticks in the plot.
            Controls the range of the x-axis display.
            Range: Depends on number of positions in the data
            Example: 5000 for 5000 positions

        xt_size: Size of x-axis tick intervals.
            Type: int
            Default: 500
            Description: Interval between x-axis tick marks.
            Controls the frequency of tick labels on the x-axis.
            Range: Typically 100 to 1000
            Example: 500 for ticks every 500 positions

        grid_width: Width of each grid position.
            Type: float
            Default: 0.01
            Description: Width of each position bar in the plot.
            Controls the visual thickness of the genotype bars.
            Range: Typically 0.005 to 0.05
            Example: 0.01 for standard bar width

    Returns:
        None. The function creates a high-resolution plot file in the
        specified format.

    Raises:
        FileNotFoundError: If genoprob_file does not exist
        ValueError: If plot parameters are invalid
        RuntimeError: If plotting fails or produces invalid output

    Notes:
        - The function automatically loads founder strain colors from
          'founder.hexcolor.info' file
        - Recombination events are counted when genotype changes between
          adjacent positions
        - The plot uses a 16x16 inch figure size for high resolution
        - Chromosomes are displayed in natural order (1, 2, ..., 19, X, Y)
        - Each chromosome shows two haplotypes with different colors
        - Recombination counts are displayed next to each chromosome
        - The plot includes a title with sample name and total recombination count
        - Output is saved at 600 DPI for publication quality
        - This function is typically run after GBRS reconstruction or interpolation
        - The visualization is essential for quality assessment and publication
    """
    if output_file is None:
        output_file = os.path.splitext(os.path.basename(genoprob_file))[0]
        output_file = f'gbrs.plotted.{output_file}.{output_format}'

    logger.info(f'Genotype Probabilities File: {genoprob_file}')
    logger.info(f'Output File: {output_file}')
    logger.info(f'Output Format: {output_format}')
    logger.info(f'Sample Name: {sample_name}')
    logger.info(f'Grid Size: {grid_size}')
    logger.info(f'XT Max: {xt_max}')
    logger.info(f'XT Size: {xt_size}')
    logger.info(f'Grid Width: {grid_width}')

    logger.info('Loading chromosome information')
    chrlens = get_chromosome_info()

    logger.info('Loading founder colors')
    hcolors = get_founder_info()
    haplotypes = hcolors.keys()
    hid = dict(zip(haplotypes, np.arange(8)))
    logger.info(f'Haplotype IDs: {hid}')

    genotypes = np.array([h1 + h2 for h1, h2 in combinations_with_replacement(haplotypes, 2)])

    #
    # Main body
    #
    logger.info(f'Loading GBRS genotype probability file: {genoprob_file}')
    genoprob = np.load(genoprob_file)

    def natural_sort(list):
        def convert(text):
            return int(text) if text.isdigit() else text.lower()

        def alphanum_key(key):
            return [convert(c) for c in re.split('([0-9]+)', key)]

        return sorted(list, key=alphanum_key)

    # intersection of ref.fa.fai chroms and those present in genoprob.
    chrs = [value for value in chrlens.keys() if value in genoprob.files]

    # natural sorting of the chrom list for display.
    chrs = natural_sort(chrs)
    num_chrs = len(chrs)

    fig = pyplot.figure()
    fig.set_size_inches((16, 16))
    ax = fig.add_subplot(111)
    ax.set_xlim(0, xt_max * grid_width + grid_width)
    ax.set_ylim(1, 95)
    num_recomb_total = 0
    for cid, c in enumerate(chrs):
        if c in genoprob.files:
            # skip drawing Y chromosome if the sample is female
            logger.debug(f'Working on {c}')
            genotype_calls = genotypes[genoprob[c].argmax(axis=0)]
            hap = []
            col1 = []
            col2 = []
            oldcol1 = 'NA'
            oldcol2 = 'NA'
            num_recomb = 0
            num_genes_in_chr = len(genotype_calls)

            for i in range(num_genes_in_chr):
                hap.append((i * grid_width, grid_width))
                c1 = hcolors[genotype_calls[i][0]]
                c2 = hcolors[genotype_calls[i][1]]

                if i > 0:
                    if c1 == c2:
                        if col1[-1] != col2[-1]:
                            # when homozygous region starts, remember the most recent het
                            oldcol1 = col1[-1]
                            oldcol2 = col2[-1]
                    else:
                        if col1[-1] == col2[-1]:
                            # when heterozygous region starts
                            if c1 == oldcol2 or c2 == oldcol1:
                                c1, c2 = c2, c1
                        elif c1 == col2[-1] or c2 == col1[-1]:
                            c1, c2 = c2, c1
                    if c1 != col1[-1] or c2 != col2[-1]:
                        num_recomb += 1

                col1.append(c1)
                col2.append(c2)
            num_recomb_total += num_recomb
            print(f'Chromsome {c} has {num_recomb} recombinations')
            # plot
            ax.broken_barh(
                hap,
                (num_chrs * 4 - cid * 4 + 1, 1),
                facecolors=col1,
                edgecolor='face',
            )
            ax.broken_barh(
                hap,
                (num_chrs * 4 - cid * 4, 1),
                facecolors=col2,
                edgecolor='face',
            )
            ax.text(
                (num_genes_in_chr + 50) * grid_width,
                num_chrs * 4 - cid * 4 + 0.5,
                f'({num_recomb})',
            )
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_right()
        ax.get_yaxis().tick_left()
        pyplot.yticks(
            ticks=np.arange(num_chrs * 4 + 1, 1, -4),
            labels=list(chrs),
            fontsize=14,
        )
        pyplot.xticks(
            ticks=np.arange(0, xt_max * grid_width, xt_size * grid_width),
            labels=[
                '%dcM' % xt
                for xt in np.arange(
                    0, xt_max * grid_size / 100, xt_size * grid_size / 100
                )
            ],
        )
        title_txt = f'Genome reconstruction: {sample_name}'
        title_txt += f'\n(Total {num_recomb_total} recombinations)'
        ax.set_title(title_txt, fontsize=18, loc='center')

    logger.info(f'Saving generated plot: {output_file}')
    fig.savefig(output_file, dpi=600, format=output_format)
    pyplot.close(fig)
    logger.info('Done')
