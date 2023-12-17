from pathlib import Path

def set_base_dir():
    """
    Set the base directory for the project.
    The function checks for the existence of the base directory in multiple locations and sets it accordingly.
    If the base directory cannot be found, a NameError is raised.
    """
    base_dir = '/Volumes/diedrichsen_data$/data/FunctionalFusion'
    if not Path(base_dir).exists():
        base_dir = '/srv/diedrichsen/data/FunctionalFusion'
    if not Path(base_dir).exists():
        base_dir = 'Y:\data\FunctionalFusion'
    if not Path(base_dir).exists():
        base_dir = '/Users/callithrix/Documents/Projects/Functional_Fusion/'
    if not Path(base_dir).exists():
        base_dir = '/Users/jdiedrichsen/Data/FunctionalFusion/'
    if not Path(base_dir).exists():
        raise (NameError('Could not find base_dir'))

    return base_dir 


def set_atlas_dir(base_dir):
    """
    Set the directory for the atlases.
    The directory is set based on the provided base directory.
    """
    atlas_dir = base_dir + f'/Atlases'
    return atlas_dir

def set_figure_dir():
    """
    Set the directory for the figures.
    The function checks for the existence of the figure directory in multiple locations and sets it accordingly.
    """
    figure_dir = "/Users/jdiedrichsen/Dropbox (Diedrichsenlab)/papers/AtlasPaper/figure_parts/"
    if not Path(figure_dir).exists():
        figure_dir = "/Users/callithrix/Dropbox/AtlasPaper/figure_parts/"

    return figure_dir

def set_export_dir(base_dir):
    """
    Set the directory for exporting.
    The directory is set based on the provided base directory.
    """
    export_dir = f'{base_dir}/../Cerebellum/ProbabilisticParcellationModel/Atlases/'
    if not Path(export_dir).exists():
        export_dir = f'{base_dir}/Atlases/'

    return export_dir

if __name__ == '__main__':
    base_dir = set_base_dir()
    atlas_dir = set_atlas_dir(base_dir)
    figure_dir = set_figure_dir()
    export_dir = set_export_dir(base_dir)