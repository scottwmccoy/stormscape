"""End-to-end example: peak-i15 over terrain for the Cameron Peak burn area.

Downloads a 10 m DEM + hillshade for a bbox around the Cameron Peak fire
(Colorado), stacks MRMS for the 2021-07-20 storm into a peak-i15 field, and
drapes the field over the hillshade. Run from anywhere once stormscape is
installed:

    python examples/example_cameron_peak.py

(or your stormscape env). Outputs land in ./example_out/.
"""

import os

from stormscape import (drape_i15, fetch_dem_and_hillshade, i15_storm_day,
                        save_fields)

OUT = "example_out"
KEY = "cameron_peak"
DATE = "2021-07-20"                       # a known post-fire storm-day
AOI = (-105.55, 40.55, -105.25, 40.80)    # W, S, E, N (lon/lat degrees)


def main():
    os.makedirs(OUT, exist_ok=True)

    # 1. DEM + hillshade (10 m 3DEP, reprojected to CONUS Albers).
    dem_path = os.path.join(OUT, f"{KEY}_dem.tif")
    hs_path = os.path.join(OUT, f"{KEY}_hillshade.tif")
    fetch_dem_and_hillshade(AOI, resolution=10, dem_path=dem_path,
                            hillshade_path=hs_path)
    print("DEM + hillshade done")

    # 2. MRMS -> peak-i15 field for the storm-day.
    res = i15_storm_day(AOI, DATE)
    save_fields(res, OUT, KEY)
    print("i15 field done:", res["meta"])

    # 3. Drape i15 over the hillshade.
    fig_path = os.path.join(OUT, f"{KEY}.png")
    drape_i15(hs_path, os.path.join(OUT, f"{KEY}_i15max.tif"),
              out_path=fig_path,
              title=f"Cameron Peak  -  peak i15 over terrain  ({DATE})")
    print("wrote", fig_path)


if __name__ == "__main__":
    main()
