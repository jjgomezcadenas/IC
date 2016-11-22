import MySQLdb

dbArt = MySQLdb.connect(host="neutrinos1.ific.uv.es", user='USER',passwd='PASS',db="NEWDB")
cursorArt = dbArt.cursor()

dbIC = MySQLdb.connect(host="neutrinos1.ific.uv.es", user='USER',passwd='PASS',db="ICNEWDB")
cursorIC = dbIC.cursor()

############
# PMT GAIN #
############
pmtGain = False
# Get latest run
cursorArt.execute('select max(MinRun) from ChannelGain where SensorID < 100')
maxRunArt = cursorArt.fetchall()[0][0]
cursorIC.execute('select max(MinRun) from PmtGain where SensorID < 100')
maxRunIC = cursorIC.fetchall()[0][0]
if maxRunArt != maxRunIC:
    pmtGain = True
    print "# PMT Gain has to be updated\n"

if pmtGain:
    # Get PMT Gains
    cursorArt.execute("select MinRun,MaxRun,SensorID,Centroid from ChannelGain where SensorID < 100 and MinRun >= {0} order by MinRun,SensorID".format(maxRunArt))
    data = cursorArt.fetchall()

    # Update old values with maxrun
    print 'UPDATE PmtGain SET MaxRun={0} where MinRun={1} and MaxRun is NULL;'.format(maxRunArt-1,maxRunIC)

    # Write new values
    for row in data:
        run = row[0]
        sensorid = row[2]
        gain = row[3]
        print 'INSERT INTO PmtGain VALUES ({0},NULL,{1},{2});'.format(run,sensorid,gain)


#############
# SIPM GAIN #
#############
sipmGain = False
# Get latest run
cursorArt.execute('select max(MinRun) from ChannelGain where SensorID > 100')
maxRunArt = cursorArt.fetchall()[0][0]
cursorIC.execute('select max(MinRun) from SipmGain where SensorID > 100')
maxRunIC = cursorIC.fetchall()[0][0]
if maxRunArt != maxRunIC:
    sipmGain = True
    print "# SiPM Gain has to be updated\n"

if sipmGain:
    # Get SiPM Gains
    cursorArt.execute("select MinRun,MaxRun,SensorID,Centroid from ChannelGain where SensorID > 100 and MinRun >= {0} order by MinRun,SensorID".format(maxRunArt))
    data = cursorArt.fetchall()

    # Update old values with maxrun
    print 'UPDATE SipmGain SET MaxRun={0} where MinRun={1} and MaxRun is NULL;'.format(maxRunArt-1,maxRunIC)

    # Write new values
    for row in data:
        run = row[0]
        sensorid = row[2]
        gain = row[3]
        print 'INSERT INTO SipmGain VALUES ({0},NULL,{1},{2});'.format(run,sensorid,gain)



#############
# SIPM MASK #
#############
# Get latest run
cursorArt.execute('select max(MinRun) from ChannelMask where SensorID > 100')
maxRunArt = cursorArt.fetchall()[0][0]
cursorIC.execute('select max(MinRun) from SipmMask where SensorID > 100')
maxRunIC = cursorIC.fetchall()[0][0]
if maxRunArt != maxRunIC:
    sipmMask = True
    print "# SiPM Mask has to be updated\n"

if sipmMask:
    # Get SiPM Gains
    #cursorArt.execute("select MinRun,MaxRun,SensorID from ChannelMask where SensorID > 100 and MinRun >= {0} order by MinRun,SensorID".format(maxRunArt))
    cursorArt.execute("select SensorID from ChannelMask where SensorID > 100 and MinRun >= {0} order by MinRun,SensorID".format(maxRunArt))
    data = cursorArt.fetchall()

    # Update old values with maxrun
    print 'UPDATE SipmMask SET MaxRun={0} where MinRun={1} and MaxRun is NULL;'.format(maxRunArt-1,maxRunIC)

    # Masked channels
    masked = []
    for row in data:
        masked.append(row[0])

    # Write new values
    for sipm in xrange(1792):
        sipmid = (sipm/64+1)*1000 + sipm % 64
        active = 0 if sipmid in masked else 1
        print 'INSERT INTO SipmMask VALUES ({0},NULL,{1},{2});'.format(maxRunArt,sipmid,active)
