# models/__init__.py

"""
========================================
Mass Estimating Models Package
========================================

This package aggregates various modules for estimating the mass of
launch vehicle engines and stages, based on several academic papers
and technical sources.

This file, __init__.py, marks the 'models' directory as a Python package.
This allows the main.py file (in the parent directory) to import modules
from this folder using syntax like:

from models import akin_mers
from models.common_params import EngineParams

Modules:
- common_params.py:      Dataclasses (EngineParams, StageParams) used by all models.
- akin_mers.py:          Akin (ENAE 791) MERs for components and SSTO analysis.
- zandbergen_engine.py:  Zandbergen (2015) statistical engine mass/size models.
- tizon_rema.py:         Tizón & Román (2017) dimensionless engine component model.
- mota_schlingloff.py:   Mota & Schlingloff (2005) analytical-statistical model.
- kibbey_stage.py:       Kibbey (2015) stage inert mass scaling model.
"""