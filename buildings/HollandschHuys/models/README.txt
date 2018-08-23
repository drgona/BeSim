Inputs (uNames):
--------------
inputs.one: 			dummy input for cost function [-].
inputs.QFlowEmbHea: 	Heat to embedded floor or ceiling [W]. Each zone has several embedded surfaces. [W]
						They can be summed up per zone. Corresponding zone-emb indices are:
						Z1: {1,8      , 15, 19}
						Z2: {2,9      , 16, 20}
						Z3: {3, 10    , 17, 21}
						Z4: {4, 11    , 18, 22}
						Z5: {23, 27   , 31, 35, 5,12, 55,57}
						Z6: {24, 28   , 32, 36, 6, 13,56, 58}
						Z7: {25, 29   , 33, 37}
						Z8: {26, 30   , 34, 38, 7, 14}
						Z9: {39, 43   , 47, 51}
						Z10:{40, 44   , 48, 52}
						Z11:{41, 45   , 49, 53}
						Z12:{42, 46  ,50, 54}
inputs.QFlowEmbCoo:		Same than inputs.QFlowEmbHea with for cooling (positive value = cooling!). [W]
						It can be summed to inputs.QFlowEmbHea using factor -1.
inputs.TSupVen[1]:		Supply ventilation temperature which is identical for all zones [K].
inputs.a:				used for cost function
inputs.slack[1]:		used for cost function
inputs.slackBou[1]:		used for cost function
winBusIn:				disturbance variables
weaBus:					disturbance variables
prescribed.TBouUp:		Upper boundary temperature for comfort. used for cost function.
prescribed.TBouLow:		Lower boundary temperature for comfort. used for cost function.
prescribed.m_flow_ven:	prescribed ventilation flow rate for each zone [kg/s].
prescribed.QGaiCon:		prescribed convective heat gain [W].
prescribed.QGaiRad:		prescribed raditiative heat gain [W].

Outputs (yNames):
--------------
All outputs starting with ineq and cost are used for the cost function and for constraints --> to be ignored.
All outputs starting with proBus1 are duplicate --> to be ignored.
TSensor:				Zone operative temperature [K]
TZoneAir:				Zone air temperature [K]
