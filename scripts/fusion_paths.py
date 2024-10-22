from pathlib import Path

def set_base_dir():
    """
    Set the base directory for the project.
    The function checks for the existence of the base directory in multiple locations and sets it accordingly.
    If the base directory cannot be found, a NameError is raised.
    """

    possible_dirs = [
        '/Volumes/diedrichsen_data$/data/',
        '/srv/diedrichsen/data',
        'Y:\data\\',
        '/Users/callithrix/Documents/Projects//',
        '/Users/jdiedrichsen/Data//',
    ]

    for directory in possible_dirs:
        if Path(directory).exists():
            return directory

    raise FileNotFoundError('Could not find base_dir')

base_dir = set_base_dir()

def set_fusion_dir(base_dir=base_dir):
    """
    Set the directory for the atlases.
    The directory is set based on the provided base directory.
    """
    fusion_dir = base_dir + f'/FunctionalFusion'
    return fusion_dir

def set_atlas_dir(base_dir=base_dir):
    """
    Set the directory for the atlases.
    The directory is set based on the provided base directory.
    """
    atlas_dir = base_dir + f'/FunctionalFusion/Atlases'
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