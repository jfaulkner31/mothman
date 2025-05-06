#########################################################################
#                               TODO BLOCK                              #
#########################################################################
# TODO: can we add descriptions to the sph / equivalence id fluxes????
# TODO change number of axial layers - crrently it doesnt pull from the mesh data but rather than user input. shouldnt ever need user input.
# TODO    can use openmc.RegularMesh() dimensuon data to fix this. until then i will be annoyed.................
# TODO SPH VALUES ON MESHES ARE COMPLETELY BUSTED - we are volume averaging in fluxes def for meshes instead of volume adding.
#########################################################################
#                  block for running input file and getting libs        #
#########################################################################
def run_input(inputFileName, doRunOpenMC, numProcessors, doMesh, doCells):
    import importlib
    import sys

    # input file name is the reference input file that we are working with.
    # doRunOpenMC is logical if we want to run openmc or not.
    # numProcessors is number of threads
    # if inputFileName[-2:] == ".py":
    #   inputFileName = inputFileName[:-2]
    inputFileName = inputFileName.replace('.py', '')
    print("now importing: "+inputFileName)
    importedFile = __import__(inputFileName) # import our input file.
    # thisModule = importlib.import_module(inputFileName, [0,0,0,0])

    # mesh based tallies
    if doMesh:
        mgxs_mesh_lib = importedFile.mgxs_mesh_lib
        mgxs_mesh_lib_transport = importedFile.mgxs_mesh_lib_transport
    else:
        mgxs_mesh_lib = 0
        mgxs_mesh_lib_transport = 0

    # cell based tallies
    if doCells:
        mgxs_cell_lib = importedFile.mgxs_cell_lib
        mgxs_cell_lib_transport = importedFile.mgxs_cell_lib_transport
    else:
        mgxs_cell_lib = 0
        mgxs_cell_lib_transport = 0

    # universe based tallies
    if doUnis:
        mgxs_uni_lib =importedFile.mgxs_uni_lib
        mgxs_uni_lib_transport = importedFile.mgxs_uni_lib_transport
    else:
        mgxs_uni_lib = 0
        mgxs_uni_lib_transport = 0

    # import groups and vols now: (i dont think we actually need vols or anything but this is here for future use.)
    fine_groups = importedFile.fine_groups
    coarse_groups = importedFile.groups
    vols = importedFile.vols # -> dictionaryt of vols from the input file for SPH factors. vols = {"f1": v1, "f2": v2, .... "g9": vg9}

    if doRunOpenMC:
        openmc.run(threads=numProcessors)
    # return some stuff that we need to work with
    return mgxs_mesh_lib, mgxs_mesh_lib_transport, mgxs_cell_lib, mgxs_cell_lib_transport, mgxs_uni_lib, mgxs_uni_lib_transport, vols, fine_groups, coarse_groups

#########################################################################
#                       get_xs block for grabbing  mesh xs              #
#########################################################################
def get_mesh_xs(mgxs_mesh_lib, mgxs_mesh_lib_transport, fine_groups, coarse_groups, homZoneArrayFile, meanType, betaLambdaType, userInputBetaArray, userInputLambdaArray, transportMethod):

    meshDim = mgxs_mesh_lib.domains[0].n_dimension # dimensions of the mesh = 3 for xyz
    meshSizing = mgxs_mesh_lib.domains[0].dimension # 27x27x2 meshwould be [27,27,2]
    if meshDim == 3:
        numAxialLayers = meshSizing[2] # if z mesh
    else:
        numAxialLayers = 1 # if xy or x mesh

    # pull pregenerated array from csv file -> see example xyArray.csv
    homZones = np.loadtxt(homZoneArrayFile, delimiter=',').astype(np.int32) # loads homZones from csv file from UI
    if meshDim == 1:
        homZones = [homZones]

    totNumberXYZones = len(homZones[0])  * len(homZones) # first dimension times 0th dimension
    totNumberZones = totNumberXYZones * numAxialLayers
    print()
    print()
    print('Total number of mesh XY Zones is: ', totNumberXYZones)
    print('Total number of mesh axial Zones is: ', numAxialLayers)
    print('Total number of mesh zones in all is: ', totNumberZones)
    print()
    print()


    #### step 1: Pull cross sections -> most expensive part...

    nuScatter = mgxs_mesh_lib.get_mgxs(mgxs_mesh_lib.domains[0], 'nu-scatter').get_condensed_xs(openmc.mgxs.EnergyGroups(coarse_groups)).get_xs() # scattering data nuScatter
    scatterProbMatrix = mgxs_mesh_lib.get_mgxs(mgxs_mesh_lib.domains[0], 'scatter probability matrix').get_condensed_xs(openmc.mgxs.EnergyGroups(coarse_groups)).get_xs()
    nuFission = mgxs_mesh_lib.get_mgxs(mgxs_mesh_lib.domains[0], 'nu-fission').get_condensed_xs(openmc.mgxs.EnergyGroups(coarse_groups)).get_xs()

    # get transport xs on our own and then condense on our own. this one has fine groups
    transport = mgxs_mesh_lib_transport.get_mgxs(mgxs_mesh_lib_transport.domains[0], 'transport').get_condensed_xs(openmc.mgxs.EnergyGroups(coarse_groups)).get_xs()
    transportCondensed = mgxs_mesh_lib_transport.get_mgxs(mgxs_mesh_lib_transport.domains[0], 'transport').get_condensed_xs(openmc.mgxs.EnergyGroups(coarse_groups)).get_xs()
    fluxFromTransport = mgxs_mesh_lib_transport.get_mgxs(mgxs_mesh_lib_transport.domains[0], 'transport').get_flux()
    # transport = mgxs_mesh_lib_transport.get_mgxs(mgxs_mesh_lib_transport.domains[0], 'transport').get_condensed_xs(openmc.mgxs.EnergyGroups(coarse_groups)).get_xs() # condense

    chiPrompt =     mgxs_mesh_lib.get_mgxs(mgxs_mesh_lib.domains[0], 'chi-prompt').get_condensed_xs(openmc.mgxs.EnergyGroups(coarse_groups)).get_xs() # chi P
    chiDelayed =    mgxs_mesh_lib.get_mgxs(mgxs_mesh_lib.domains[0], 'chi-delayed').get_condensed_xs(openmc.mgxs.EnergyGroups(coarse_groups)).get_xs() # chi delayed
    beta =          mgxs_mesh_lib.get_mgxs(mgxs_mesh_lib.domains[0], 'beta').get_condensed_xs(openmc.mgxs.EnergyGroups(coarse_groups)).get_xs() # beta
    decayRate =     mgxs_mesh_lib.get_mgxs(mgxs_mesh_lib.domains[0], 'decay-rate').get_condensed_xs(openmc.mgxs.EnergyGroups(coarse_groups)).get_xs()
    velocityInverse =    mgxs_mesh_lib.get_mgxs(mgxs_mesh_lib.domains[0], 'inverse-velocity').get_condensed_xs(openmc.mgxs.EnergyGroups(coarse_groups)).get_xs()
    redAbs =        mgxs_mesh_lib.get_mgxs(mgxs_mesh_lib.domains[0], 'reduced absorption').get_condensed_xs(openmc.mgxs.EnergyGroups(coarse_groups)).get_xs()
    flux =          mgxs_mesh_lib.get_mgxs(mgxs_mesh_lib.domains[0], 'reduced absorption').get_condensed_xs(openmc.mgxs.EnergyGroups(coarse_groups)).get_flux()
    kappa_fission = mgxs_mesh_lib.get_mgxs(mgxs_mesh_lib.domains[0], 'kappa-fission').get_condensed_xs(openmc.mgxs.EnergyGroups(coarse_groups)).get_xs()

    #### step 2: we make a dictionary with keys indicating the zones from the mask above.

    thisDict = {}
    for axialLayerNumber in range(numAxialLayers):
        for yLine in homZones[::-1]:
            for xVal in yLine:
                k = "axialLayer_"+str(axialLayerNumber+1)+"_cellZone_"+str(xVal) # starting index at 0 but turn into 1 in naming conventions
                thisDict[k] = []

    #### step 3: append iterators to the array (iterators are indexes of the array)

    iterator = 0
    # now we append iterators
    for axialLayerNumber in range(numAxialLayers):
        for yLine in homZones[::-1]:  # iterate backwards -> ensures that we start at lower left corner and move in +x and +y direction.
            for xVal in yLine:
                k = "axialLayer_"+str(axialLayerNumber+1)+"_cellZone_"+str(xVal)
                thisDict[k].append(iterator)
                iterator += 1
    # we now have a dictionary with keys and indexes where keys are the
    # HOMOGENOUS ZONES and INDX is WHERE IN THE XS ARRAY THOSE ZONES ARE

    # now we want to go key by key and average out everything.
      # now set keys for each reaction type...

    nuFissionDict = {}
    transportDict = {}
    removalDict = {}
    diffusionDict = {}
    chiPromptDict = {}
    chiDelayedDict = {}
    betaDict = {}
    decayRateDict = {}
    velocityDict = {}
    nuSTimesProbDict = {}
    fluxDict = {}
    kappa_fissionDict = {}

    #### step 4: make a loop that goes from lower left corner to top right
        ## start at yMin zMin xMin
        ## go in +x then +y until all of this layer are done then +z
    for key in thisDict.keys():
        lengthForAvg = len(thisDict[key])
        print()
        print()
        print("Now creating / homogenizing XS dictionary for mesh region -", key, "with number of homogenous regions =", lengthForAvg)
        # initialize cross sections as arrays...
        nuFissionArr     = []
        transportArr     = []
        removalArr       = []
        diffusionArr     = []
        chiPromptArr     = []
        chiDelayedArr    = []
        betaArr          = []
        decayRateArr     = []
        velocityArr      = []
        nuSTimesProbArr  = []
        fluxArr          = []
        kappa_fissionArr = []
        fluxArrFromTransport = []

        for xsIndex in thisDict[key]:
            # NOW WE GRAB XS FROM TALLIES AND PLACE THOSE XS INTO A BIG LIST
            # scattering matrix  - hardest one to do...
            nuSTimesProb = [] # this turns into scattering matrix

            diagOfScatterProbScatter_XS = []
            nuScatter_XS         = nuScatter[xsIndex] # nu*scatter xs
            scatterProbMatrix_XS = scatterProbMatrix[xsIndex] # scattering prob matrix

            index = 0
            for row in scatterProbMatrix_XS:
                nuSTimesProb.append(row*nuScatter_XS[index])
                diagOfScatterProbScatter_XS.append(row[index]*nuScatter_XS[index])
                index+=1

            nuSTimesProb = np.vstack(nuSTimesProb) # new edition - vstack turns list of arrays into a stacked array FINAL SCTTR MTRX

            nuSTimesProbArr.append(nuSTimesProb)
            transportArr.append(transport[xsIndex])
            diffusionArr.append( 1/3/transport[xsIndex] )
            nuFissionArr.append(     nuFission[xsIndex]             )
            chiPromptArr.append(     chiPrompt[xsIndex]             )
            chiDelayedArr.append(    chiDelayed[xsIndex]            )
            betaArr.append(          beta[xsIndex]                  )
            decayRateArr.append(     decayRate[xsIndex]             )
            velocityArr.append(      1/velocityInverse[xsIndex]     )
            fluxArr.append(          flux[xsIndex]                  )
            kappa_fissionArr.append( kappa_fission[xsIndex]         )
            fluxArrFromTransport.append(fluxFromTransport[xsIndex]            )
            removalArr.append(  redAbs[xsIndex] + nuScatter_XS - diagOfScatterProbScatter_XS    )

        # now once we get a list of all xs in the index, we can simply take the mean of them
        # TODO check that we need mean of fluxes? we might actually need sum of we are doing equivalence procedure........
        # put output values into dicts  -> output values are averages amongst the lists that we have made ....
        if meanType == 0:
            nuFissionDict[key]       = np.mean(nuFissionArr, axis=0)
            raise Exception("error due to transport mean not being done correctly!! need to change the meanType option and this error before user can use this option! always use meanType = 1 for now")
            transportDict[key]     = np.mean(transportArr, axis=0)

            removalDict[key]       = np.mean(removalArr, axis=0)
            diffusionDict[key]     = np.mean(diffusionArr, axis=0)
            chiPromptDict[key]     = np.mean(chiPromptArr, axis=0)
            chiDelayedDict[key]    = np.mean(chiDelayedArr, axis=0)
            velocityDict[key]      = np.mean(velocityArr, axis=0)
            nuSTimesProbDict[key] = np.mean(nuSTimesProbArr, axis=0)
            fluxDict[key]          = np.mean(fluxArr, axis=0)
            kappa_fissionDict[key] = np.mean(kappa_fissionArr, axis=0)
        elif meanType == 1:
            nuSTimesProbDict[key]  = fluxWeightedAverageXSMatrix(fluxArr, nuSTimesProbArr) # ok
            nuFissionDict[key]     = fluxWeightedAverageXS(fluxArr, nuFissionArr) # ok

            if transportMethod == 0:
                transportDict[key]     = fluxWeightedTransportXS(fluxArrFromTransport, transportArr, coarse_groups, fine_groups) # ok
            elif transportMethod == 1:
                transportDict[key] = transportCondensed[xsIndex]

            removalDict[key]       = fluxWeightedAverageXS(fluxArr, removalArr)   # ok
            diffusionDict[key]     = inverseWeightedAverageXS(fluxArr, diffusionArr) # ok ## TODO fix errr here shouldnt be inverse weighted????
            velocityDict[key]      = inverseWeightedAverageXS(fluxArr, velocityArr) #ok
            fluxDict[key]          = fluxAverage(fluxArr) # ok
            kappa_fissionDict[key] = fluxWeightedAverageXS(fluxArr, kappa_fissionArr) # ok
            chiPromptDict[key]     = chiPromptAverageXS(fluxArr, chiPromptArr, nuFissionArr) # ok
            chiPromptDummyPlugin = chiPromptDict[key] # dummy variable to pass into the chiDelayed processing - see error in chidelayed for description of why we do this.
            chiDelayedDict[key]    = chiDelayedAverageXS(fluxArr, betaArr, nuFissionArr, chiDelayedArr, chiPromptDummyPlugin, key) # ok but watch out for nan errors in sparse regions
        else:
            raise Exception("Cross section mean averaging type not known (meanType=?)")
        if betaLambdaType == 0:
            betaDict[key]          = np.mean(betaArr, axis=0)
            decayRateDict[key]     = np.mean(decayRateArr, axis=0)
        elif betaLambdaType == 1:
            betaDict[key]          = betaAverageXS(fluxArr, betaArr, nuFissionArr, key) # ok
            decayRateDict[key]     = decayRateAverageXS(fluxArr, betaArr, nuFissionArr, decayRateArr) #ok
        elif betaLambdaType == 2:
            betaDict[key]          = userInputBetaArray
            decayRateDict[key]     = userInputLambdaArray
        else:
            raise Exception("Beta Lambda Type (betaLambdaType) not known.")

    # now make reaction rate csv files
    makeRxRateFile(meshDim, meshSizing, kappa_fission, flux, "kappaFissionFlux") # makes a reaction rate file over the tallied mesh for output

    return nuFissionDict, transportDict, removalDict, diffusionDict, chiPromptDict, chiDelayedDict, betaDict, decayRateDict, velocityDict, nuSTimesProbDict, fluxDict, kappa_fissionDict

#########################################################################
#                       flux weighted averaging                         #
#########################################################################
def fluxWeightedAverageXSMatrix(fluxArr, xsArr):
    numGroups = len(fluxArr[0]) # returns number of groups [1,2] woudl return 2
    outputArr = np.zeros([numGroups, numGroups])
    for thisGroup in range(numGroups): # 0, 1, 2, ... numGroups - 1
        fluxSummed = 0.0
        xsRow = np.zeros(numGroups) # creates empty array of all zeros that is the size of numGroups
        for zoneNum, regionFlux in enumerate(fluxArr): # region flux is the flux array([....]) of this region
            # print("xsRow is:", xsRow)
            # print("and adding", xsArr[zoneNum][thisGroup], "times", fluxArr[zoneNum][thisGroup], "=",xsArr[zoneNum][thisGroup] * fluxArr[zoneNum][thisGroup])
            xsRow += xsArr[zoneNum][thisGroup] * fluxArr[zoneNum][thisGroup] # rx rate for this zone and this group added to the xsRow
            fluxSummed += fluxArr[zoneNum][thisGroup] # summed the flux for this zone and this group

        outputArr[thisGroup] = xsRow/fluxSummed # append this homogenized row to the scattering xs array

    return outputArr  # xsArrOutput

def fluxWeightedTransportXS(fluxArr, transportArr, coarse_groups, passed_fine_groups):
    coarse_groups = coarse_groups[::-1] # reverse groups since the input needs them in forward direction (low to high) and out needs high to low
    passed_fine_groups = passed_fine_groups[::-1]
    coarse_groups = np.array(coarse_groups)
    #### first condense all the xs from fine groups to coarse groups
    indexes = np.array([])
    for group in coarse_groups:
        # print(group)
        indexes = np.append(indexes, np.argwhere(passed_fine_groups == group))

    # coarse transport xs and fluxes
    coarse_transport = np.zeros(len(coarse_groups)-1)
    coarse_flux_sums = np.zeros(len(coarse_groups)-1)
    low = 0 ## always start at index of 0 - > usually 2e7 value
    for i, fineIndex in enumerate(indexes[1::]):
    # i is the index of the indexes , fineIndex is the index in the fine_structure e.g. 0,1,3,5,251....
        fineIndex = int(fineIndex) # convert the index into an int for indexing.
        for regionIndex, regionXS in enumerate(transportArr):
            regionFlux = fluxArr[regionIndex]
            multiply = regionFlux[low:fineIndex] / regionXS[low:fineIndex] #### divide since we are doing transport XS (not multiply) - same as condensing diffusion coeffs.
            fluxesToAdd = regionFlux[low:fineIndex]
            for m, ratio in enumerate(multiply): # if nan when we divide tr (e.g. transport xs has value of 0)
                if math.isnan(ratio):
                    multiply[m] = 0.0
                    # fluxesToAdd[m] = 0.0
                if multiply[m] < 0.0: # if somehow negative - set to 0 - doesnt contribute to the tally.
                    multiply[m] = 0.0
                    # fluxesToAdd[m] = 0.0


            coarse_transport[i] = coarse_transport[i] + np.sum(multiply)
            coarse_flux_sums[i] = coarse_flux_sums[i] + np.sum(fluxesToAdd)
            # print(low, fineIndex, np.sum(multiply))

        low = fineIndex # reset the fine index so next time around we go g-1 -> g

    coarse_transport = coarse_transport / coarse_flux_sums
    coarse_transport=  1 / coarse_transport # divide again at the end... since we did inverse weighting and we want to return the coasre_transport




    return coarse_transport

def fluxWeightedAverageXS(fluxArr, xs):
    numGroups = len(fluxArr[0])
    output = np.zeros([numGroups])
    for thisGroup in range(numGroups):
        fluxSummed = 0.0
        rxRate = 0.0
        for zoneNum, regionFlux in enumerate(fluxArr):
            rxRate += xs[zoneNum][thisGroup] * fluxArr[zoneNum][thisGroup]
            fluxSummed += fluxArr[zoneNum][thisGroup]

        output[thisGroup] = rxRate / fluxSummed

    return output

def chiPromptAverageXS(fluxArr, xs, weightArr):
    # weight arr is a weight to be also multiplied with the fluxes.
    # for example, chi_g = sum(chi_g_vi * nsf_g * phi_g) / sum(nsf_g * phi_g)
    # in this example, nsf_g would be the new weightArr to be used in the weighting.

    numGroups = len(fluxArr[0])
    output = np.zeros([numGroups])
    for thisGroup in range(numGroups):
        weightedSummed = 0.0
        rxRate = 0.0
        for zoneNum, regionFlux in enumerate(fluxArr):
            rxRate += xs[zoneNum][thisGroup] * fluxArr[zoneNum][thisGroup] * weightArr[zoneNum][thisGroup]
            weightedSummed += fluxArr[zoneNum][thisGroup] * weightArr[zoneNum][thisGroup]

        output[thisGroup] = rxRate / weightedSummed

    return output

def inverseWeightedAverageXS(fluxArr, xs):
    # uses the inverse of the xs as the weight - returns the non-inverse e.g. uses 1/v * phi but returns: v

    numGroups = len(fluxArr[0])
    output = np.zeros([numGroups])
    for thisGroup in range(numGroups):
        weightedSummed = 0.0
        rxRate = 0.0
        for zoneNum, regionFlux in enumerate(fluxArr):
            rxRate += 1/xs[zoneNum][thisGroup] * fluxArr[zoneNum][thisGroup]
            weightedSummed += fluxArr[zoneNum][thisGroup]

        output[thisGroup] = 1 / (rxRate / weightedSummed) # vBAR = (1/v * phi) / sum(phi)

    return output

def betaAverageXS(fluxArr, xs, weightArr, key):
    #
    numGroups = len(fluxArr[0])
    numDelayedGroups = len(xs[0])
    output = np.zeros([numDelayedGroups])

    for d in range(numDelayedGroups): # d is the DNP index

        # print(d)
        beta_d = 0.0
        for zoneNum, beta_arr in enumerate(xs):
            # comppute total nsF * phi * beta_d_g over all groups.
            for g in range(numGroups):
                beta_d += fluxArr[zoneNum][g] * weightArr[zoneNum][g] * xs[zoneNum][d][g]

        output[d] = beta_d

    rxnRate = 0
    for zoneNum, thisFlux in enumerate(fluxArr):
        for g in range(numGroups):
            rxnRate += fluxArr[zoneNum][g] * weightArr[zoneNum][g]

    for d in range(numDelayedGroups):
        output[d] = output[d] / rxnRate

    if np.isnan(output).any():
        print("!!!!!!!!!!!!!!!!!!!!!!!!! NAN ERROR NAN ERROR NAN ERROR !!!!!!!!!!!!!!!!!!!!!!!!!")
        print("NAN error thrown for region / key =", key, "- this means that a NAN was detected in this region. This is due to this region having total FISSION rx rate of 0 for all homogenized zones and thus triggering a division by 0.")
        print("Note that this error may actually be valid for some regions where there are no fissions.")
        print("changing all beta values to 0 now.")
        print("!!!!!!!!!!!!!!!!!!!!!!!!! NAN ERROR NAN ERROR NAN ERROR !!!!!!!!!!!!!!!!!!!!!!!!!")
        output = np.zeros(numDelayedGroups)

    return output

def chiDelayedAverageXS(fluxArr, betaXS, nsfXS, chiXS, chiPromptXS, key):
    numGroups = len(fluxArr[0])
    numDelayedGroups = len(betaXS[0])
    output = np.zeros([numDelayedGroups, numGroups]) # d rows and g cols
    numZones = len(fluxArr)
    # do top part:
    for d in range(numDelayedGroups):
        topSum = np.zeros(numGroups)
        for g in range(numGroups):
            groupSum = 0.0
            for zoneNum in range(numZones):
                groupSum += betaXS[zoneNum][d][g] * nsfXS[zoneNum][g] * fluxArr[zoneNum][g] * chiXS[zoneNum][d][g]

            topSum[g] = groupSum

        output[d] = topSum

    # do bottom part:
    for d in range(numDelayedGroups):
        rxRate = np.zeros(numGroups)
        for g in range(numGroups):
            for zoneNum in range(numZones):
                rxRate[g] = rxRate[g] + betaXS[zoneNum][d][g] * nsfXS[zoneNum][g] * fluxArr[zoneNum][g] * sum(chiXS[zoneNum][d]) # sum() will be either 1 or 0. if 0 - then add no weight

        output[d] = output[d] / rxRate
        output[d] = output[d] / sum(output[d]) # normalize just in case - double checked and we really dont need to but normalizing helps numerically anywasy
    if np.isnan(output).any():
        print("!!!!!!!!!!!!!!!!!!!!!!!!! NAN ERROR NAN ERROR NAN ERROR !!!!!!!!!!!!!!!!!!!!!!!!!")
        print("NAN error thrown for region / key =", key, "- this means that a NAN was detected in this region. This is due to this region having total FISSION rx rate of 0 for all homogenized zones and thus triggering a division by 0.")
        print("The delayed neutron spectrum for this region will now be changed to the PROMPT spectrum for this region!!!")
        print("This means that all user defined dnp spectra will be the same as the prompt spectra (chi_d_g = chi_g for all values of d)")
        print("Note that this error will be OK if it occurs in non-fissile regions!!!")
        print("if this error is occuring in fissile regions then it is likely due to undersampling and poor statistics for the delayed neutrons and should be validated.")
        print("!!!!!!!!!!!!!!!!!!!!!!!!! NAN ERROR NAN ERROR NAN ERROR !!!!!!!!!!!!!!!!!!!!!!!!!")
        # raise Exception("nan error detected in chidelayed arrays")
        for idx,row in enumerate(output):
            output[idx] = chiPromptXS

    return output

def decayRateAverageXS(fluxArr, betaXS, nsfXS, lambdaXS):
    numGroups = len(fluxArr[0])
    numDelayedGroups = len(betaXS[0])
    numZones = len(fluxArr)
    output = np.zeros(numDelayedGroups)
    for d in range(numDelayedGroups):
        summedTop = 0.0
        summedBot = 0.0
        for zoneNum in range(numZones):
            for g in range(numGroups):
                summedTop += lambdaXS[zoneNum][d] * nsfXS[zoneNum][g] * fluxArr[zoneNum][g] * betaXS[zoneNum][d][g]
                summedBot += nsfXS[zoneNum][g] * fluxArr[zoneNum][g] * betaXS[zoneNum][d][g]

        output[d] = summedTop / summedBot

    ### throw error here if any of the decay consts = 0

    return output

def fluxAverage(fluxArr):
    # returns the groupwise summed fluxes across all macro regions

    return np.sum(fluxArr, axis=0)

#########################################################################
#                 flux weighted rx rates files                          #
#########################################################################
def makeRxRateFile(meshDim, meshSizing, rxRateInput, flux, rxRateNameForFile):
    # makes an output csv file using the rxRateInput (a cross section usually) and the flux.
    # meshes are arranged so the FIRST ARRAY in the csv file (e.g. the one appearing at the very top) corresponds to the lowest elevation set of voxels.
    # otherwise, the xy plane reads normally
    thisValIndex = 0
    if meshDim == 3:
        # 3dimensional mesh
        numX = meshSizing[0]
        numY = meshSizing[1]
        numZ = meshSizing[2]
        arr = np.zeros([numZ,numY,numX]) # numZ arrays of numY rows and numX cols
        for thisZ in range(numZ):
            for thisY in range(numY):
                for thisX in range(numX):
                        # values but starting at lower left and working our way up.
                        # if numY meshes is 27 then numY - 1 - thisY will be 26 to start - e.g. the very bottom of the array.
                        # doing this for Y and Z guarantees we move rightwards in X, up in Y, and up in Z
                        arr[numZ - 1 - thisZ, numY - 1 - thisY, thisX] = np.sum(rxRateInput[thisValIndex] * flux[thisValIndex]) # fill array with corresponding data
                        thisValIndex += 1
        dummyArray = ['']*numX
        emptyArr = np.zeros([numX]) # empty arrayof zeros
        for thisZ in arr:
            emptyArr = np.vstack((emptyArr, thisZ))
            emptyArr = np.vstack((emptyArr, dummyArray))
        emptyArr = emptyArr[1::] # get rid of first dummy array in z direction
        csvName = rxRateNameForFile+"_"+str(numX)+"x"+str(numY)+"x"+str(numZ)+".csv"
        pd.DataFrame(emptyArr).to_csv(csvName, index=False) # then make a csv!!!!

    elif meshDim == 2:
        numX = meshSizing[0]
        numY = meshSizing[1]
        numZ = 0
        arr = np.zeros([numY, numX])
        for thisY in range(numY):
            for thisX in range(numX):
                # values but starting at lower left and working our way up.
                arr[numY - 1 - thisY, thisX] = np.sum(rxRateInput[thisValIndex] * flux[thisValIndex]) # fill array with corresponding data
                thisValIndex += 1
        csvName = rxRateNameForFile+"_"+str(numX)+"x"+str(numY)+"x"+str(numZ)+".csv"
        pd.DataFrame(emptyArr).to_csv(csvName, index=False) # then make a csv!!!!

    elif meshDim == 1: # 1-D mesh
        numX = meshSizing[0]
        numY = 0
        numZ = 0
        arr = np.zeros([numX])
        for thisX in range(numX):
            # values but starting at lower left and working our way up.
            arr[thisX] = np.sum(rxRateInput[thisValIndex] * flux[thisValIndex]) # fill array with corresponding data
            thisValIndex += 1
        csvName = rxRateNameForFile+"_"+str(numX)+"x"+str(numY)+"x"+str(numZ)+".csv"
        pd.DataFrame(arr).to_csv(csvName, index=False) # then make a csv!!!!

    return

#########################################################################
#                       get_xs block for grabbing cell xs               #
#########################################################################
def get_cell_xs(this_mgxs_lib, mgxs_lib_transport, domainNumber, fine_groups, coarse_groups, betaLambdaType, userInputBetaArray, userInputLambdaArray, domainName, transportMethod):
    i = domainNumber
    mgxs_libs = [this_mgxs_lib]

    nuScatter = mgxs_libs[0].get_mgxs(mgxs_libs[0].domains[i], 'nu-scatter').get_condensed_xs(openmc.mgxs.EnergyGroups(coarse_groups)).get_xs() # scattering data nuScatter
    scatterProbMatrix = mgxs_libs[0].get_mgxs(mgxs_libs[0].domains[i], 'scatter probability matrix').get_condensed_xs(openmc.mgxs.EnergyGroups(coarse_groups)).get_xs()
    nuSTimesProb = []
    diagOfScatterProbScatter = [] # diagonal of scater prob matrix * nuScatter -> used for removal XS
    index = 0
    for row in scatterProbMatrix:
        nuSTimesProb.append(row*nuScatter[index])
        diagOfScatterProbScatter.append(row[index]*nuScatter[index])

        index+=1

    # now we have to get fission data and transport and delayed data.
    nuFission = mgxs_libs[0].get_mgxs(mgxs_libs[0].domains[i], 'nu-fission').get_condensed_xs(openmc.mgxs.EnergyGroups(coarse_groups)).get_xs() # fission data

    # print("tally num ----- ", i)
    transportPre = mgxs_lib_transport.get_mgxs(mgxs_lib_transport.domains[i], 'transport').get_condensed_xs(openmc.mgxs.EnergyGroups(coarse_groups)).get_xs() #### .get_condensed_xs(openmc.mgxs.EnergyGroups(coarse_groups)) # transport data
    transportCondensed = mgxs_lib_transport.get_mgxs(mgxs_lib_transport.domains[i], 'transport').get_condensed_xs(openmc.mgxs.EnergyGroups(coarse_groups)).get_xs()

    transportFlux = mgxs_lib_transport.get_mgxs(mgxs_lib_transport.domains[i], 'transport').get_flux() ####
    if transportMethod == 0:
        transport = fluxWeightedTransportXS([transportFlux], [transportPre], coarse_groups, fine_groups)
    elif transportMethod == 1:
        transport = transportCondensed

    diffusion = 1/3/transport
    chiPrompt = mgxs_libs[0].get_mgxs(mgxs_libs[0].domains[i], 'chi-prompt').get_condensed_xs(openmc.mgxs.EnergyGroups(coarse_groups)).get_xs() # chi P
    chiDelayed = mgxs_libs[0].get_mgxs(mgxs_libs[0].domains[i], 'chi-delayed').get_condensed_xs(openmc.mgxs.EnergyGroups(coarse_groups)).get_xs() # chi delayed
    velocity = 1/mgxs_libs[0].get_mgxs(mgxs_libs[0].domains[i], 'inverse-velocity').get_condensed_xs(openmc.mgxs.EnergyGroups(coarse_groups)).get_xs()
    redAbs = mgxs_libs[0].get_mgxs(mgxs_libs[0].domains[i], 'reduced absorption').get_condensed_xs(openmc.mgxs.EnergyGroups(coarse_groups)).get_xs()
    flux = mgxs_libs[0].get_mgxs(mgxs_libs[0].domains[i], 'reduced absorption').get_condensed_xs(openmc.mgxs.EnergyGroups(coarse_groups)).get_flux()
    kappa_fission = mgxs_libs[0].get_mgxs(mgxs_libs[0].domains[i], 'kappa-fission').get_condensed_xs(openmc.mgxs.EnergyGroups(coarse_groups)).get_xs()

    removal = redAbs + nuScatter - diagOfScatterProbScatter

    one_group = openmc.mgxs.EnergyGroups()
    one_group.group_edges = np.array([coarse_groups[0], coarse_groups[-1]])
    if betaLambdaType == 0 | betaLambdaType == 1:

        # beta = mgxs_libs[0].get_mgxs(mgxs_libs[0].domains[i], 'beta').get_condensed_xs(one_group).get_xs() # new beta with I entries
        beta = betaAverageXS(np.array([flux]),    np.array([ mgxs_libs[0].get_mgxs(mgxs_libs[0].domains[i], 'beta').get_xs() ]),   np.array([nuFission])  , domainName)
        decayRate = mgxs_libs[0].get_mgxs(mgxs_libs[0].domains[i], 'decay-rate').get_xs()

    elif betaLambdaType == 2:
        beta = userInputBetaArray
        decayRate = userInputLambdaArray
    else:
        raise Exception("beta lambda type unknown from user input? (betaLambdaType=???)")
    #beta = mgxs_libs[0].get_mgxs(mgxs_libs[0].domains[i], 'beta').get_xs() # beta - old beta with I*G entries


    return nuFission, transport, removal, diffusion, chiPrompt, chiDelayed, beta, decayRate, velocity, nuSTimesProb, flux, kappa_fission


#########################################################################
#                       get_xs block for grabbing univ xs               #
#########################################################################
def get_uni_xs(this_mgxs_lib, mgxs_lib_transport, domainNumber, fine_groups, coarse_groups, betaLambdaType, userInputBetaArray, userInputLambdaArray, domainName, transportMethod):
    i = domainNumber
    mgxs_libs = [this_mgxs_lib]

    nuScatter = mgxs_libs[0].get_mgxs(mgxs_libs[0].domains[i], 'nu-scatter').get_condensed_xs(openmc.mgxs.EnergyGroups(coarse_groups)).get_xs() # scattering data nuScatter
    scatterProbMatrix = mgxs_libs[0].get_mgxs(mgxs_libs[0].domains[i], 'scatter probability matrix').get_condensed_xs(openmc.mgxs.EnergyGroups(coarse_groups)).get_xs()
    nuSTimesProb = []
    diagOfScatterProbScatter = [] # diagonal of scater prob matrix * nuScatter -> used for removal XS
    index = 0
    for row in scatterProbMatrix:
        nuSTimesProb.append(row*nuScatter[index])
        diagOfScatterProbScatter.append(row[index]*nuScatter[index])

        index+=1

    # now we have to get fission data and transport and delayed data.
    nuFission = mgxs_libs[0].get_mgxs(mgxs_libs[0].domains[i], 'nu-fission').get_condensed_xs(openmc.mgxs.EnergyGroups(coarse_groups)).get_xs() # fission data

    # print("tally num ----- ", i)
    transportPre = mgxs_lib_transport.get_mgxs(mgxs_lib_transport.domains[i], 'transport').get_condensed_xs(openmc.mgxs.EnergyGroups(coarse_groups)).get_xs() #### .get_condensed_xs(openmc.mgxs.EnergyGroups(coarse_groups)) # transport data
    transportFlux = mgxs_lib_transport.get_mgxs(mgxs_lib_transport.domains[i], 'transport').get_flux() ####
    transportCondensed = mgxs_lib_transport.get_mgxs(mgxs_lib_transport.domains[i], 'transport').get_condensed_xs(openmc.mgxs.EnergyGroups(coarse_groups)).get_xs()

    if transportMethod == 0:
        transport = fluxWeightedTransportXS([transportFlux], [transportPre], coarse_groups, fine_groups)
    elif transportMethod == 1:
        transport = transportCondensed

    diffusion = 1/3/transport
    chiPrompt = mgxs_libs[0].get_mgxs(mgxs_libs[0].domains[i], 'chi-prompt').get_condensed_xs(openmc.mgxs.EnergyGroups(coarse_groups)).get_xs() # chi P
    chiDelayed = mgxs_libs[0].get_mgxs(mgxs_libs[0].domains[i], 'chi-delayed').get_condensed_xs(openmc.mgxs.EnergyGroups(coarse_groups)).get_xs() # chi delayed
    decayRate = mgxs_libs[0].get_mgxs(mgxs_libs[0].domains[i], 'decay-rate').get_condensed_xs(openmc.mgxs.EnergyGroups(coarse_groups)).get_xs()
    velocity = 1/mgxs_libs[0].get_mgxs(mgxs_libs[0].domains[i], 'inverse-velocity').get_condensed_xs(openmc.mgxs.EnergyGroups(coarse_groups)).get_xs()
    redAbs = mgxs_libs[0].get_mgxs(mgxs_libs[0].domains[i], 'reduced absorption').get_condensed_xs(openmc.mgxs.EnergyGroups(coarse_groups)).get_xs()
    flux = mgxs_libs[0].get_mgxs(mgxs_libs[0].domains[i], 'reduced absorption').get_condensed_xs(openmc.mgxs.EnergyGroups(coarse_groups)).get_flux()
    kappa_fission = mgxs_libs[0].get_mgxs(mgxs_libs[0].domains[i], 'kappa-fission').get_condensed_xs(openmc.mgxs.EnergyGroups(coarse_groups)).get_xs()

    removal = redAbs + nuScatter - diagOfScatterProbScatter

    one_group = openmc.mgxs.EnergyGroups()
    one_group.group_edges = np.array([coarse_groups[0], coarse_groups[-1]])
    if betaLambdaType == 0 | betaLambdaType == 1:

        # beta = mgxs_libs[0].get_mgxs(mgxs_libs[0].domains[i], 'beta').get_condensed_xs(one_group).get_xs() # new beta with I entries
        beta = betaAverageXS(np.array([flux]),    np.array([ mgxs_libs[0].get_mgxs(mgxs_libs[0].domains[i], 'beta').get_xs() ]),   np.array([nuFission])  , domainName)
        decayRate = mgxs_libs[0].get_mgxs(mgxs_libs[0].domains[i], 'decay-rate').get_xs()

    elif betaLambdaType == 2:
        beta = userInputBetaArray
        decayRate = userInputLambdaArray
    else:
        raise Exception("beta lambda type unknown from user input? (betaLambdaType=???)")
    #beta = mgxs_libs[0].get_mgxs(mgxs_libs[0].domains[i], 'beta').get_xs() # beta - old beta with I*G entries


    return nuFission, transport, removal, diffusion, chiPrompt, chiDelayed, beta, decayRate, velocity, nuSTimesProb, flux, kappa_fission


#########################################################################
#                      Import and grabbing other stuff                  #
#########################################################################

## decalre this in function definition - make expandable in future. add anything that is hardcoded to this block for later on when we generalize.
import copy
import sys
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import openmc
import copy
import os
from   collections import OrderedDict
import h5py
from   itertools import product
import numpy as np
from   xml.etree import ElementTree as ET
from   getfilename import *
import openmc
from openmc._xml import clean_indentation
import openmc
import os
import h5py
import matplotlib.pyplot as plt
import openmc.checkvalue as cv
import os
import h5py
import pandas as pd
import math
from openmc._xml import clean_indentation
from openmcToGriffinMakeTallies import *

#########################################################################
#                       User input settings                             #
#########################################################################
# user inputs for the user to change: numGenerationsInInput, doRunOpenMC, numProcessesToRunOpenMC,
#                                     LibraryName, NGroup, num_delayed_groups, Generator, TimeCreated,
#                                     sphOn, doDiffusion, outputXMLName, temperature
numGenerationsInInput = 20      # number in statepoint e.g. statepoint.numGenerationsInInput.h5 -> e.g. statepoint.500.h5 we would put 500 for this variable
doRunOpenMC = False              # turn on to true if we want to run the openmc file as well - almost always false
numProcessesToRunOpenMC = 15     # num processes to run openmc with if we do...
numStatesToTabulate = 1 # number of tabulation states
LibraryName = 'MSRE_XS_LIB_1D_CORE' # title of library
NGroup = 12 # number groups
num_delayed_groups = 6
Generator = "jmf" # name of person making file
TimeCreated = "11-4-2023-1441est" # time created.
sphOn = True                    # generates equivalence data for SPH fluxes.
 # this is if there is multiple calcs and therefore multiple libraries to choose from.
tabulation = OrderedDict()
tabulation['Tfuel'] = [922,]
tabulation['Tgraphite'] = [922,]
legendre_order = 0 # use L=1 for diffusion calculations since that is what all the XML's on VTB use # however griffin errors when L=1 so use L=0 for now ig.
doDiffusion = False #### if false then DOES NOT write diffusion coeffs to xml file
outputXMLName = "TESTING_ONLY.xml" # name of output xml
doMesh = True    # true if there are mesh tallies to use
doCells = True   # true if there are cell tallies to use
doUnis =  False   # true if there are universe tallies to use in xs gen
inputFileName = "msre_core_1d.py" # name of openmc input files. DO NOT INCLUDE .py EXTENSION!!!!! . also make sure to comment out the openmc.run in the input file.
fileNewName = "statepoint."+str(numGenerationsInInput)+".h5"     # name of statepoint we want
homZoneArrayFile = 'array_1d.csv'            # array file for mesh xs tallies - we can leave empty if we arent doing mesh tallies. see xyArray.csv for an example
condenseFromFine = 1        # new option in r3 - if == 1 then we condense ALL xs from fine energy grids.
transportMethod = 1         # 1 = openmc default weighting. 0 = custom weighted inverse.
meanType = 1                # type of mean to take when averaging mesh cells - 0 for simple average - 1 for flux weighted
betaLambdaType = 1          # type of mean to take when averaging beta nad lambda - 0 for simple average - 1 for flux weighted (see meeting minutes) and 2 for user defined
userInputBetaArray = []     # when defining xs, mandates that we use a specific array of beta values in the xs process if we are using betaLambdaType = 2 for user defined values
userInputLambdaArray = []   # when defining xs, mandates that we use a specific array of lambda vlaues in the xs process if we are using betaLambdaType = 2 for user defined values
axialLayerMultiplier = 1000 # used for mesh tally processing. if we have cell label (from the xyArray.csv file) of 15 and axial layer number 8 then the ID will be
                            # 8 * axialLayerMultiple + 15. so if axialLayMult = 1000 -> we get 8015 as the ID of the xs with name of axialLayer_8_cellZone_15.
                            # this guy should be larger than the element number in the xy array. so if we have >1000 elements in xy array then we would need to use 10000 as the multiplier

print("File name = ", fileNewName)
#########################################################################
#########################################################################
#########################################################################

mgxs_mesh_lib, mgxs_mesh_lib_transport, mgxs_cell_lib, mgxs_cell_lib_transport , mgxs_uni_lib, mgxs_uni_lib_transport, vols, fine_groups, coarse_groups = run_input(inputFileName, doRunOpenMC, numProcessesToRunOpenMC, doMesh, doCells)

NGroup = len(coarse_groups) - 1

print()
print()
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("~~~~~~~~~~~~~~~~~~ Input file successfully ran! ~~~~~~~~~~~~~~~~~~~~~")
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

print("~~~~~~~~~~~~~~~~~ Now running xs generation tool! ~~~~~~~~~~~~~~~~~~~")
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print()
print("Now loading statepoint file:", fileNewName)
SP = openmc.StatePoint(fileNewName) # open statepoint for the filename

if (doMesh):
    print("\n\n\n\n\nUser indicated doMesh to be true meaning there are mesh tallies to write - ")
    print("Now loading mesh normal xs library ...")
    mgxs_mesh_lib.load_from_statepoint(SP) # loads mesh tallies for general XS
    print("Now loading mesh transport xs library ...")
    mgxs_mesh_lib_transport.load_from_statepoint(SP) # loads mesh transport cross section tallies

if (doCells):
    print("User indicated doCells to be true meaning there are cell tallies to write - ")
    print("Now loading cell based transport xs library ...")
    mgxs_cell_lib_transport.load_from_statepoint(SP)
    print("Now loading cell based normal xs library ...")
    mgxs_cell_lib.load_from_statepoint(SP)
if (doUnis):
    print("User indicated doUnis to be True - meaning that there are universe tallies to write")
    print("Now loading universe based transport xs library ...")
    mgxs_uni_lib_transport.load_from_statepoint(SP)
    print("Now loading universe based normal xs library ...")
    mgxs_uni_lib.load_from_statepoint(SP)

#############################################
######### start of file #####################
#############################################
xml = ET.Element("OpenMC_2_Griffin")
root = ET.SubElement(xml, "Multigroup_Cross_Section_Libraries")
root.set("Name", str(LibraryName))
root.set("NGroup", str(NGroup))
if sphOn:
    equivalence = ET.SubElement(xml, "Equivalence_Data_Library")
    equivalence.set("Name", LibraryName)
    equivalence.set("NGroup", str(NGroup))


##############################################
######### Sort through domains now  ##########
#  we will have combo of cell / mesh domains #
##############################################
# print(mgxs_mesh_lib_transport.domains[0]) # was here to print out RegularMesh
domain_number = -1
cell_domain_index = -1 # set to -1 and only add once we are done sorting through mesh domains
uni_domain_index = -1

if doMesh:
    nuFissionDict,transportDict,removalDict,diffusionDict,chiPromptDict,chiDelayedDict, betaDict,decayRateDict,velocityDict,nuSTimesProbDict,fluxDict, kappa_fissionDict = get_mesh_xs(mgxs_mesh_lib, mgxs_mesh_lib_transport,
                                                                                                                                                                                       fine_groups, coarse_groups,
                                                                                                                                                                                       homZoneArrayFile,
                                                                                                                                                                                       meanType, betaLambdaType,
                                                                                                                                                                                       userInputBetaArray, userInputLambdaArray,
                                                                                                                                                                                       transportMethod)
    num_mesh_cells = len(nuFissionDict.keys())
    meshKeyList = list(nuFissionDict.keys())
else:
    num_mesh_cells = 0
if doCells:
    num_cell_domains = np.size(mgxs_cell_lib.domains[:])
else:
    num_cell_domains = 0
if doUnis:
    num_uni_domains = np.size(mgxs_uni_lib.domains[:])
else:
    num_uni_domains = 0

print("\nTotal number of cell domains = ", num_cell_domains)
print("Total number of mesh voxel domains = ", num_mesh_cells)
print("Total number of universe domains = ", num_uni_domains)

prevID = 0
while domain_number < -1 + num_mesh_cells + num_cell_domains + num_uni_domains:  ### get size of mesh - e.g. number of mesh cells plus number of cell domains
    # print( num_mesh_cells)
    domain_number += 1 # iterable for domain number that we iterate through
    prevID += 1 # adds 1 to prevID -> if we are past the mesh writing stage, this just indexes +1 for every material ID that we iterate through.
                # so if meshes end at 20110, then the next ID will be 20111
                # if we do not do any meshing, then this very naturally will have ID=1 for the first cell/universe that we go through!!!
    ID = prevID #
    # doMesh condition
    if (domain_number < num_mesh_cells) & doMesh:
        name = meshKeyList[domain_number] # grabs domain from mesh key list
        print('Now doing mesh number: '+str(ID), "with name:", name)
        nuFission = nuFissionDict[name]
        transport = transportDict[name]
        removal = removalDict[name]
        diffusion = diffusionDict[name]
        chiPrompt = chiPromptDict[name]
        chiDelayed = chiDelayedDict[name]
        beta = betaDict[name]
        decayRate = decayRateDict[name]
        velocity = velocityDict[name]
        scattering = nuSTimesProbDict[name]
        flux = fluxDict[name]
        kappa_fission = kappa_fissionDict[name]

        axialLayerNumber = int(name[11:name.find('_', 11, len(name))]) ##### k = "axialLayer_"+str(axialLayerNumber)+"_cellZone_"+str(xVal)
        cellZoneNumber =   int( name[   (name.find('_', 11, len(name)) + len('_cellZone_')) :len(name)   ] ) ####
        ID = int (  (axialLayerNumber) * axialLayerMultiplier + cellZoneNumber  )# override ID -> add 1 to the axialLayerNumber index for convenience
        prevID = max(prevID, ID) # override prev id if we are still doing meshes....
                                # need to do max since meshes are not necessarily sorted in numerical order. this creates a bug since prevID used further down will overwrite mesh_ids
                                # need to guarantee that prevID will ALWAYS be the max(meshID) + 1 at the end of printing meshes.

    # doCells condition
    elif (domain_number >= num_mesh_cells) & (domain_number < num_cell_domains + num_mesh_cells) & doCells:
        cell_domain_index += 1
        thisDomain = mgxs_cell_lib.domains[cell_domain_index]
        name = thisDomain.name
        nuFission, transport, removal, diffusion, chiPrompt, chiDelayed, beta, decayRate, velocity, scattering, flux, kappa_fission = get_cell_xs(mgxs_cell_lib, mgxs_cell_lib_transport,
                                                                                                                                                  cell_domain_index, fine_groups,
                                                                                                                                                  coarse_groups,
                                                                                                                                                  betaLambdaType,userInputBetaArray, userInputLambdaArray, name,
                                                                                                                                                  transportMethod)
        print('Now doing cell domain with ID = ', str(ID), " and cell domain name = ", name)

    # doUniverse conditions
    elif (domain_number >= num_cell_domains + num_mesh_cells) & doUnis:
        uni_domain_index += 1
        thisUni = mgxs_uni_lib.domains[uni_domain_index]
        name = thisUni.name
        nuFission, transport, removal, diffusion, chiPrompt, chiDelayed, beta, decayRate, velocity, scattering, flux, kappa_fission = get_uni_xs(mgxs_uni_lib, mgxs_uni_lib_transport,
                                                                                                                                                  uni_domain_index, fine_groups, coarse_groups,
                                                                                                                                                  betaLambdaType, userInputBetaArray, userInputLambdaArray, name,
                                                                                                                                                  transportMethod)
        print("Now doing universe with ID = ", str(ID), "and universe domain name = ", name)


    library = ET.SubElement(root, "Multigroup_Cross_Section_Library")
    library.set("ID", str(ID)) # sets ID -> this is what is used in ISOXML / GRIFFIN to map XS at the end of the day
    library.set("Description", name) # this just provides user description
    library.set("Ver", str(1.0))
    library.set("Generator", Generator)
    library.set("TimeCreated", TimeCreated)
    GridCoordNames = ET.SubElement(library, "Tabulation") # tabulation vairable.
    ReferenceGridIndex = ET.SubElement(library, "ReferenceGridIndex")


    if sphOn:
      equivalence_data = ET.SubElement(equivalence, "EquivalenceData")
      equivalence_data.set("ID", str(ID))

    GridCoordNames.text = ""
    ReferenceGridIndex.text = ""
    IndexDict = OrderedDict()

    if sphOn:
        GridCoordNames2 = ET.SubElement(equivalence_data, "Tabulation")
        ReferenceGridIndex2 = ET.SubElement(equivalence_data, "ReferenceGridIndex")

    for key in tabulation:
        IndexDict[key] = {}
        ReferenceGridIndex.text = ReferenceGridIndex.text+str(len(tabulation[key]))+" "

        GridCoordNames.text = GridCoordNames.text+key+" "
        variable = ET.SubElement(library, key)
        if sphOn:
            variable2 = ET.SubElement(equivalence_data, key)
        values = ""
        countvalue = 0
        for value in tabulation[key]:
            countvalue = countvalue+1
            IndexDict[key].update({value: countvalue})
            values = values+str(value)+" "
        variable.text = values[:-1]
        if sphOn:
            variable2.text = values[:-1]
    ReferenceGridIndex.text = ReferenceGridIndex.text[:-1]
    GridCoordNames.text = GridCoordNames.text[:-1]

    if sphOn:
      GridCoordNames2.text = GridCoordNames.text
      ReferenceGridIndex2.text = ReferenceGridIndex.text


    ## all reaction types:
    rnxs = ['chi-prompt', 'chi-delayed','nu-fission',
            'transport','nu-scatter matrix','chi-delayed','beta','decay-rate']

    ignoreFissionStuff = all(a == 0 for a in nuFission) # returns true if the entire nuFission array is values of zero - necessary for new update of isoxml


    allReactions = ET.SubElement(library, "AllReactions")
    tablewiseReactions = ET.SubElement(library, "TablewiseReactions")

    if ignoreFissionStuff: # truncate severely if we want to ignore all the fission stuff- needed as of November 2023 for griffin due to recoverable energy changes?
                           # notr eally sure exact date when change was but code now errors due to kappa_fission being 0 when there is a fissionXS tabulated (even if it is 0)
        allReactions.text = "Removal NeutronVelocity Transport Scattering"
        tablewiseReactions.text = "Removal NeutronVelocity Transport Scattering"
    else:
        allReactions.text = "Removal DNFraction FissionSpectrum DNSpectrum DNPlambda NeutronVelocity nuFission kappaFission Transport Scattering"
        tablewiseReactions.text = "Removal DNFraction FissionSpectrum DNSpectrum DNPlambda NeutronVelocity nuFission kappaFission Transport Scattering"

    for i, f in enumerate([fileNewName]): # i is index f is filename from fn.listfn
        Table = ET.SubElement(library, "Table")
        Table.set("gridIndex", ReferenceGridIndex.text)
        if sphOn:
            Table2 = ET.SubElement(equivalence_data, "Table")
            Table2.set("gridIndex", ReferenceGridIndex.text)
            Flux = ET.SubElement(Table2, "Flux")
            Flux.text = ""
            # SP = openmc.StatePoint(fn.listfn[i])
            # t = SP.get_tally(name=name) # use name=name for proper sph tally. originally was error here but no longer? hopefully....
            # thisTally = t.mean[::-1] # get tally and reverse it.
            #thisTally = t.mean[:] # nonreversed.
            for result in flux: ### this is where sph factors are written out!
                Flux.text = Flux.text+str(result)+" " # use /vols[name] if we weant to normalize to vol....
                # print(vols[name])
                # print(str(result))
            Flux.text = Flux.text[:-1]
        Isotope = ET.SubElement(Table, "Isotope")
        Isotope.set("Name", "pseudo")
        Isotope.set("L", str(legendre_order))
        Isotope.set("I", str(num_delayed_groups))

        Tablewise = ET.SubElement(Table, "Tablewise")
        Tablewise.set("L", str(legendre_order))
        Tablewise.set("I", str(num_delayed_groups))
        datafile = datafile = h5py.File(f,'r')
        #### now go through reactions one by one and print out to XML
        # Removal DNFraction FissionSpectrum DNSpectrum DNPlambda NeutronVelocity nuFission Transport Scattering


        #### REMOVAL ######
        Reaction = ET.SubElement(Tablewise, "Removal")
        Reaction.set('index', 'g')
        Reaction.text = ""
        for d in removal:
            Reaction.text = Reaction.text+str(d)+" "
        ### TRANSPORT #####
        Reaction = ET.SubElement(Tablewise, "Transport")
        Reaction.set('index', 'g')
        Reaction.text = ""
        for d in transport:
            Reaction.text = Reaction.text+str(d)+" "


        if not(ignoreFissionStuff):
            ### nu-Fission #####
            Reaction = ET.SubElement(Tablewise, "nuFission")
            Reaction.set('index', 'p')
            Reaction.text = ""
            for d in nuFission:
                Reaction.text = Reaction.text+str(d)+" "

            ### kappa-fission ###
            Reaction = ET.SubElement(Tablewise, "kappaFission")
            Reaction.set('index', 'p')
            Reaction.text = ""
            for d in kappa_fission:
                Reaction.text = Reaction.text+str(d)+" "

        ### Diffusion #####
        if doDiffusion:
            Reaction = ET.SubElement(Tablewise, "Diffusion")
            Reaction.set('index', 'g')
            Reaction.text = ""
            for d in diffusion:
                Reaction.text = Reaction.text+str(d)+" "

        ### velocity #####
        Reaction = ET.SubElement(Tablewise, "NeutronVelocity")
        Reaction.set('index', 'g')
        Reaction.text = ""
        for d in velocity:
            Reaction.text = Reaction.text+str(d)+" "


        if not(ignoreFissionStuff):
            ### chi-prompt #####
            Reaction = ET.SubElement(Tablewise, "FissionSpectrum")
            Reaction.set('index', 'g')
            Reaction.text = ""
            for d in chiPrompt:
                Reaction.text = Reaction.text+str(d)+" "

            ### chi-delayed #####
            Reaction = ET.SubElement(Tablewise, "DNSpectrum")
            Reaction.set('index', 'gi')
            Reaction.text = ""
            for d in chiDelayed:
                for d2 in d:
                    Reaction.text = Reaction.text+str(d2)+" "

            ### DNFRACTIONS #####
            Reaction = ET.SubElement(Tablewise, "DNFraction")
            Reaction.set('index', 'i')
            Reaction.text = ""
            for d2 in beta:
                Reaction.text = Reaction.text+str(d2)+" "

            ### DNPLAMBDA #####
            Reaction = ET.SubElement(Tablewise, "DNPlambda")
            Reaction.set('index', "i")
            Reaction.text = ""
            for d in decayRate:
                Reaction.text = Reaction.text+str(d)+" "

        #### SCATTERING MATRIX ########
        scattering = np.transpose(scattering)
        Reaction = ET.SubElement(Tablewise, "Scattering")
        Reaction.set("index", "pgl")
        Reaction.set("profile", "1")
        Profile  = ET.SubElement(Reaction, "Profile")
        Profile.text = ""
        ng = str(NGroup)
        for indxer in range(NGroup):
            Profile.text = Profile.text+"\n\t\t\t\t\t\t\t"+"1 "+ng
        Profile.text = Profile.text+"\n\t\t\t\t\t\t"
        Value  = ET.SubElement(Reaction, "Value")
        Value.text = ""
        for row in scattering:
            Value.text = Value.text + "\n\t\t\t\t\t\t\t"
            for val in row:
                myString = "{:<25}".format(str(val))
                Value.text = Value.text + myString
        Value.text = Value.text + "\n\t\t\t\t\t\t"





##############################
##### cleanup and write ######
##############################
print("Now Making XML File:")
clean_indentation(xml)
tree = ET.ElementTree(xml)
tree.write(outputXMLName, xml_declaration=True, encoding='utf-8')
##############################

################################
