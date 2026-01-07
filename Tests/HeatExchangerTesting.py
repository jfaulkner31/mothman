import HeatExchanger.HX as HX

area = 1.0
htc = 1000.0
temp_rel_tol = 1e-10

hot_dict = {
  'pressure': 101325.0,
  'mdot': 1.0
}

cold_dict = {
  'pressure': 101325.0,
  'mdot': 1.0
}


the_hx = HX.HeatExchangerObject(area=area, htc=1000.0, hot_dict=hot_dict, cold_dict=cold_dict)
