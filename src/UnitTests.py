if __name__ == "__main__":
    # These imports are used despite what ruff thinks
    from Bloop.PythoniseMathematica import PythoniseMathematicaUnitTests # noqa: F401
    from Bloop.ParsedExpression import ParsedExpressionUnitTests # noqa: F401
    from Bloop.TransitionFinder import TransitionFinderUnitTests # noqa: F401
    from Bloop.Z2_ThreeHiggsBmGenerator import BmGeneratorUnitTests # noqa: F401
    from Bloop.PDGData import PDGUnitTests # noqa: F401

    from unittest import main

    main()
