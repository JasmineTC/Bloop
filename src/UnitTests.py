if __name__ == "__main__":
    # These imports are used despite what ruff thinks
    from ThreeHiggs.PythoniseMathematica import PythoniseMathematicaUnitTests # noqa: F401
    from ThreeHiggs.ParsedExpression import ParsedExpressionUnitTests # noqa: F401
    from ThreeHiggs.TransitionFinder import TransitionFinderUnitTests # noqa: F401
    from ThreeHiggs.Z2_ThreeHiggsBmGenerator import BmGeneratorUnitTests # noqa: F401
    from ThreeHiggs.PDGData import PDGUnitTests # noqa: F401

    from unittest import main

    main()
