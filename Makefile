PRINT	   = pr

RM	      = /bin/rm -f

TARGET	  = photon

EXTRAOBJS       = 

##############################################

all:	    $(TARGET)       

install:	$(TARGET)

$(TARGET):      $(TARGET).o
	$(CXX) $(LDFLAGS) $(TARGET).o -o $@ $(LIBS)

$(TARGET).o : $(TARGET).c $(TARGET).h
	$(CC) -DMAIN $(CFLAGS) -c $(TARGET).c

baremetal: photon.c photon.h
	$(CC) -DMAIN $(CFLAGS) $(TARGET).c -ffreestanding -nostdlib -static 

clean:  
	$(RM) $(TARGET).o $(TARGET)
##############################################
