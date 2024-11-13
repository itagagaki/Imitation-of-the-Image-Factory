ALL :
	cd Compiler ; $(MAKE)
	cd Renderer ; $(MAKE)

clean :
	cd Compiler ; $(MAKE) clean
	cd Renderer ; $(MAKE) clean
