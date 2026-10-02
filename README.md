# PDAF Binding for TIE-GCM

This repository contains the coupling layer between the
[Parallel Data Assimilation Framework (PDAF)](https://github.com/PDAF/PDAF)
and the [NCAR HAO](https://ncar.ucar.edu/) [Thermosphere Ionosphere Electrodynamics General Circulation Model (TIE-GCM)](https://github.com/NCAR/tiegcm).
It is a dependency of the
[**TIE-GCM-PDAF fork**](https://github.com/rainbowsend/TIE-GCM-PDAF), which
compiles it together with TIE-GCM into a single executable. It contains no
makefiles and cannot be built or run on its own; to perform TIE-GCM
simulations with PDAF, refer to the fork.

The source code was written by Armin Corbin at the
Astronomical, Physical and Mathematical Geodesy Group,
[Institute for Geodesy and Geoinformation](https://www.igg.uni-bonn.de/apmg/de),
[University of Bonn](https://www.uni-bonn.de/en/home?set_language=en).

## License

> [!IMPORTANT]
> **This software is part of the NCAR TIE-GCM. Use is governed by the Open
> Source Academic Research License Agreement contained in the file
> [`tiegcmlicense.txt`](tiegcmlicense.txt)**
> (also [online](https://github.com/NCAR/tiegcm?tab=License-1-ov-file)).
>
> That agreement permits use for **research, academic, and non-profit
> purposes only**. In particular:
>
> - **No commercial use.** Neither the software, nor any work containing it,
>   nor data products generated with it may be sold, licensed, or transferred
>   for a fee (§3a).
> - **No operational use.** It may not be used in applications intended to
>   create forecast, nowcast, or hindcast products (§3a).
> - **Derivative works** must carry change notices and be distributed on an
>   open-source basis (§3b).

The source code in this repository is licensed under the
[GNU Lesser General Public License v3.0](LICENSES/LGPL-3.0-or-later.txt), matching the
license of [PDAF](https://github.com/PDAF/PDAF?tab=LGPL-3.0-1-ov-file).
Because the code can only be used together with TIE-GCM, the TIE-GCM license
applies in addition; the file [`LICENSE`](LICENSE) summarizes how the two
relate.

This repository contains neither TIE-GCM nor PDAF source code and grants no
rights to use, modify, or redistribute either; users must obtain them
separately under their own terms. This project is not affiliated with,
endorsed by, or owned by NCAR. Some files (the PDAF-side callback routines,
e.g. `init_pdaf.F90`, `collect_state_pdaf.F90`, `distribute_state_pdaf.F90`)
are adapted from PDAF's tutorial and template code.
