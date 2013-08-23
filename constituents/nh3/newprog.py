import ta
evenNumbers = []
for i in range(100):
    evenNumbers.append(i*2.)

oddNumbers = []
for i in evenNumbers:
    oddNumbers.append(i+1)

whatIsEvenAndOddAddedTogether = []
for i,e in enumerate(evenNumbers):
    whatIsEvenAndOddAddedTogether.append(evenNumbers[i]+oddNumbers[i])

