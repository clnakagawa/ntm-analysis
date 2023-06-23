import getTest, mafftTest, queryCheck, time, os

def main():
    for i in range(20):
        getTest.main()
        cwd = os.getcwd()
        mafftTest.run()
        time.sleep(60)
        os.chdir(cwd)
        queryCheck.main()

if __name__ == "__main__":
    main()