art_make( BASENAME_ONLY
    LIBRARY_NAME  LArHSN_larhsn_ScanRecoSelectionParameters
    MODULE_LIBRARIES

    LArHSN_larhsn_ScanRecoSelectionParameters

		larcore_Geometry
	  larcore_Geometry_Geometry_service
	  larcoreobj_SummaryData

	  larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
	  larsim_MCCheater_BackTracker_service

		lardata_Utilities
		lardata_DetectorInfoServices_DetectorPropertiesServiceStandard_service
		lardataobj_RawData
	  lardataobj_RecoBase
	  lardataobj_AnalysisBase
	  lardata_RecoObjects larreco_Deprecated
		larreco_RecoAlg

		larevt_Filters
		nusimdata_SimulationBase
	 	larpandora_LArPandoraInterface
	  
	  ${LARRECO_LIB}
	  ${LARDATA_LIB}

	  ${ART_FRAMEWORK_CORE}
	  ${ART_FRAMEWORK_PRINCIPAL}
	  ${ART_FRAMEWORK_SERVICES_REGISTRY}
	  ${ART_FRAMEWORK_SERVICES_OPTIONAL}
	  ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}

	  art_Persistency_Common
	  art_Persistency_Provenance
	  art_Utilities

	  ${MF_MESSAGELOGGER}
	  ${MF_UTILITIES}

	  ${FHICLCPP}
	  ${CETLIB}

	  ${ROOT_GEOM}
	  ${ROOT_XMLIO}
	  ${ROOT_GDML}
	  ${ROOT_BASIC_LIB_LIST}
		${ROOT_CORE}
		${ROOT_EG}
    )

install_headers()
install_fhicl()
install_source()
