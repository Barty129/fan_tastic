import fantastic.turbomachinery as turbo


def main():
    inlet = turbo.Node(101325, 298, 287)
    outlet = turbo.Node(101325, 298, 287)

    compressor = turbo.Compressor()
    turbine = turbo.Turbine()
    vacuum = turbo.Vacuum()

    rig = turbo.Rig(inlet, [compressor, turbine, vacuum], outlet)
    rig.add_shaft([0, 1])



if __name__ == "__main__":
    main()
