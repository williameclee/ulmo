import copernicusmarine

copernicusmarine.subset(
    dataset_id="cmems_mod_glo_phy_my_0.083deg_P1M-m",
    variables=["so"],
    minimum_longitude=-180,
    maximum_longitude=179.9166717529297,
    minimum_latitude=-80,
    maximum_latitude=90,
    start_datetime="1993-01-01T00:00:00",
    end_datetime="2025-11-01T00:00:00",
    minimum_depth=0.49402499198913574,
    maximum_depth=5727.9169921875,
)
