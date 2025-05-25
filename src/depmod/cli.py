from __future__ import annotations

import click

import depmod
from depmod.log import logger

@click.command()
def invoke_depmod_cli():
    logger.info(f"DEformation Path for Molecular Dynamics v{depmod.__version__}")
