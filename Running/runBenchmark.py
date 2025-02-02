from ThreeHiggs.UserInput import UserInput
userinput = UserInput()
args = userinput.parse()

from ThreeHiggs.UserInput import Stages
if args.firstStage <= Stages.convertMathematica <= args.lastStage:
    from ThreeHiggs.convertMathematica import convertMathematica
    convertMathematica(args)

if args.firstStage <= Stages.minimization <= args.lastStage:
    from ThreeHiggs.minimization import minimization
    minimization(args)

