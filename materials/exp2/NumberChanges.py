from psychopy import visual, core, data, event, logging, gui, misc
import numpy, random, os

if not os.path.exists('data/'):
    os.makedirs('data/')

# settings
stimSize = .75
locBound = [9.8, 7.3]
locSep = 2
SameKey = 'm'
ChangeKey = 'z'
quitKey = 'q'
frameRate = 60
RT = core.Clock()
Ns = [5,8]
trialsPerBlock = 68

# 4 change conditions (1:4 changes)
# 2 set sizes (randomly intermixed)
# 34 trials per change/same, N, and N_changes = 8 blocks of 68 trials = 544 overall

# stimuli
colours = [[-1,-1,-1], # 'black', 
          [ 1, 1, 1], # 'white', 
          [ 1,-1,-1], # 'red', 
          [-1,-1, 1], # 'blue', 
          [-1, 0,-1], # 'green', 
          [ 1, 1,-1], # 'yellow', 
          [ 1,.3,-1], # 'orange', 
          [-1, 1, 1], # 'cyan', 
          [ 0,-1, 0], # 'purple', 
          [-1, 0, 0], # 'darkbluegreen' / teal...
          [ 1,-1, 1], # fuchsia
          [-1, 1,-1], # lime
          [ 0,-1,-1]] # maroon
Nstim = len(colours)

# most stimuli needed for the N = 8, N_c = 4 condition. 
# should aim for 13 colours... (fuchsia, lime, brown or maroon)

widthLim = [-locBound[0]/2, locBound[0]/2]
heightLim = [-locBound[1]/2, locBound[1]/2]

# this function handles the main trial sequence
def trialSeq(block):
    nDone = 0
    for thisTrial in block:
        # # # SET-UP TRIAL (done before time critical stuff)
        # extract relevant information for trial
        N = thisTrial['SetSize']
        ttype = thisTrial['SoC']
        nChange = thisTrial['nChange']
        
        # create trial stimuli and masks
        trialMasks = []
        trialStim = []
        for s in range(N):
            trialStim.append(visual.Rect(win, width = stimSize, height = stimSize))
            trialMasks.append(visual.GratingStim(win, tex=maskTex, mask=None, size=stimSize, sf=1.0/stimSize, texRes=256))
        
        # select trial locations
        trialLocs = []
        for i in range(N):
            works = 0
            while works < 1: # while we havent got a location that works (i.e. is not more than 2 deg from others...)
                candidate = [random.uniform(widthLim[0], widthLim[1]), random.uniform(heightLim[0], heightLim[1])] # generate random candidate
                close = 0 # test is this close
                if numpy.sqrt((candidate[0])**2 + (candidate[1])**2) > locSep: # if candidate is far from centre
                    for loc in trialLocs:
                       dist = numpy.sqrt((loc[0] - candidate[0])**2 + (loc[1] - candidate[1])**2)
                       if dist < locSep: # if the candidate location is near an already sampled one...
                           close += 1
                else:
                    close += 1
                if close == 0: # if the candidate is not close to any
                   trialLocs.append(candidate) # add to trial list
                   works += 1 # this one works...
        
        # place trial stimuli in location
        for s in range(N):
            trialStim[s].setPos(trialLocs[s])
            trialMasks[s].setPos(trialLocs[s])
        
        # select study stimuli (indices to store later)
        study = random.sample(range(Nstim), N)
        
        changeItems = random.sample(range(N), nChange)
        newItems = random.sample([x for x in range(Nstim) if x not in study], nChange)
        
        test = study[:] # copy
        if ttype == 'c': # if this is a change trial
            for (i, new) in zip(changeItems, newItems):
                test[i] = new
        
        # # # MAIN TRIAL SEQUENCE
        if nChange == 1:
            wd = 'change'
        else:
            wd = 'changes'
        instrText.setText('Ready - %s possible %s' %(nChange, wd))
        
        # 1000 ms warning (with reminder of change probability)
        for fr in range(frameRate):
            instrText.draw()
            win.flip()
        
        # 1000 ms fixation
        for fr in range(frameRate):
            fix.draw()
            win.flip()
        
        # set study colours
        for s in range(N):
            trialStim[s].setLineColor(colours[study[s]])
            trialStim[s].setFillColor(colours[study[s]])
        
        # present study array for 500 ms
        for fr in range(frameRate/2):
            [trialStim[i].draw() for i in range(N)]
            win.flip()
        
        # 500 ms blank screen
        for fr in range(frameRate/2):
            win.flip()
        
        # present masks for 500 ms
        for fr in range(frameRate/2):
            [trialMasks[i].draw() for i in range(N)]
            win.flip()
        
        # set test colours + draw
        for s in range(N):
            trialStim[s].setLineColor(colours[test[s]])
            trialStim[s].setFillColor(colours[test[s]])
        [trialStim[i].draw() for i in range(N)]
        
        RT.reset() # start RT counter and present test
        win.flip()
        
        # # # COLLECT RESPONSE
        event.clearEvents(eventType='keyboard')
        keyP = 0
        while keyP == 0:
            for key in event.getKeys():
                if key in [SameKey]:
                    ans = 0
                    keyP += 1
                elif key in [ChangeKey]:
                    ans = 1
                    keyP += 1
                elif key in [quitKey]:
                    core.quit()
        # to test - comment above lines and uncomment below
         
        #core.wait(.5); ans = random.choice([1, 0])
        
        reactionTime = RT.getTime()
        
        nDone += 1
        
        # # # WRITE TO CSV
        dataFile.write('%i,%i,%.4f,%i,%i,%.4f,%s,%s\n' %(nDone, N, nChange, int(ttype == 'c'), ans, reactionTime, '-'.join(str(x) for x in study), '-'.join(str(x) for x in test)))


# this function generates a block of trials
def genBlock(SetSizes, nChange, nTrials):
    trialsPerSS = nTrials/len(SetSizes)
    if trialsPerSS % 2 != 0: # if there's not an even number of trials
        core.quit()
    trialList = []
    for ss in SetSizes:
        for st in range(trialsPerSS/2):
            trialList.append({'SetSize': ss, 'SoC': 's', 'nChange': nChange})
        for dt in range(trialsPerSS/2):
            trialList.append({'SetSize': ss, 'SoC': 'c', 'nChange': nChange})
    random.shuffle(trialList) # randomise order
    return trialList


def instr(message):
    event.clearEvents()
    t = visual.TextStim(win, text=message, color = "black", height = 1, pos = [0,-4], wrapWidth = 25)
    t.draw()
    win.flip()
    keypress = 0
    while keypress < 1:
        for key in event.getKeys():
            if key in [quitKey]:
                core.quit()
            else:
                keypress += 1

def main():
    # open up data file + gui
    expInfo = {'Participant' : 0, 'Gender' :['Male','Female'], 'Age' : 0}
    expInfo['dateStr'] = data.getDateStr()
    
    dlg = gui.DlgFromDict(expInfo, title = "Participant Info", fixed = ['dateStr'], order=['Participant', 'Age', 'Gender'])
    if dlg.OK:
        LogFile = "NC_participantInfo"
        infoFile = open('data/' + LogFile+'.csv', 'a')
        infoFile.write('%s,%s,%s,%s\n' %(expInfo['Participant'], expInfo['Gender'], expInfo['Age'], expInfo['dateStr']))
        infoFile.close()
    else:
        core.quit()
    
    fileName = 'Participant' + str(expInfo['Participant']) + '_NC' # + expInfo['dateStr']
    global dataFile, win, patchCols, maskTex, fix, instrText, instrImage
    dataFile = open('data/' + fileName+'.csv', 'w')
    dataFile.write('Trial, SetSize, Cond, testChange, respChange, RT, study, test\n')
    
    win = visual.Window([1024, 768], fullscr=True, allowGUI= False, units = 'deg', monitor = 'testMonitor')
    patchCols = [colours[col] for col in [4, 2, 8, 6, 12, 7, 11, 1, 9, 0, 5, 10, 3, 2, 5, 1]] # random.sample([colours[col] for col in range(Nstim)*2], 16)
    maskTex = numpy.array(patchCols).reshape((4,4,3))

    fix = visual.TextStim(win, text = '+', color = "black", height = 1, pos = [0,0])
    instrText = visual.TextStim(win, color = "black", height = 1, pos = [0,0], wrapWidth = 30) # before each trial
    instrImage = visual.ImageStim(win, contrast = 1, pos = [0,5]) # before each block
    
    blockList = []
    for i in range(2):
        for j in random.sample([1, 2, 3, 4], 4):
            blockList.append(j)
    
    # main experimental trials
    for block in range(len(blockList)):
        nc = blockList[block]
        trials = genBlock(Ns, blockList[block], trialsPerBlock)
        # instruction screen + image 
        instrImage.setImage('images/changes' + str(nc) + '.png')
        instrImage.draw()
        if nc == 1:
            sqText = 'square'
        else:
            sqText = 'squares'
        instr("In this block if a change occurs %s %s will change colour between study and test. There will be an equal number of change and same trials.\n\n\nPress %s for SAME or %s for CHANGE\n\n\nPress any key to start." %(nc, sqText, SameKey, ChangeKey))
        win.flip(); core.wait(.5)
        trialSeq(trials)
        if block < len(blockList) - 1:
            instr('End of block\n\n\nPress any key to begin next block')
        else:
            instr('End of Experiment\n\n\nThank you for taking part!')

# run
main()

