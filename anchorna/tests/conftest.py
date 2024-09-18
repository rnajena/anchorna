import pytest


def pytest_addoption(parser):
    parser.addoption('--full', action='store_true',
                     default=False, help='run the full test set including slow tests')

def pytest_collection_modifyitems(config, items):
    if not config.getoption('--full'):
        skip_slow = pytest.mark.skip(reason='slow test: need --full option to run')
        for item in items:
            if 'slowtest' in item.keywords:
                item.add_marker(skip_slow)
    # explicitly add filter warnings to markers so that they have a higher
    # priority than command line options, e.g. -W error
    for item in items:
        for fwarn in config.getini('filterwarnings'):
            item.add_marker(pytest.mark.filterwarnings(fwarn))
