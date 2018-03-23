plot(mod::PhyloNetworks.TraitSubstitutionModel) = error("plot not defined for $(typeof(mod)).")

"""
    plot(model::TwoBinaryTraitSubstitutionModel)

Display substitution rates for a trait evolution model for two
possibly dependent binary traits; using `R` via `RCall`.
Adapted from fitPagel functions found in the `R` package `phytools`.

# Examples

```
julia> using PhyloNetworks

julia> m = TwoBinaryTraitSubstitutionModel([2.0,1.2,1.1,2.2,1.0,3.1,2.0,1.1],
                ["carnivory", "noncarnivory", "wet", "dry"])
Substitution model for 2 binary traits, with rate matrix:
                     carnivory-wet    carnivory-dry noncarnivory-wet noncarnivory-dry
    carnivory-wet                *           1.0000           2.0000           0.0000
    carnivory-dry           3.1000                *           0.0000           1.1000
 noncarnivory-wet           1.2000           0.0000                *           2.0000
 noncarnivory-dry           0.0000           2.2000           1.1000                *

julia> plot(m);
```
"""
function plot(object::PhyloNetworks.TwoBinaryTraitSubstitutionModel)
    R"""
    signif<-3
    plot.new()
    par(mar=c(1.1,2.1,3.1,2.1))
    plot.window(xlim=c(0,2),ylim=c(0,1),asp=1)
    """
    R"""
    mtext("Two Binary Trait Substitution Model",side=3,adj=0,line=1.2,cex=1.2)
    arrows(x0=0.15,y0=0.15,y1=0.85,lwd=2,length=0.1)
    arrows(x0=0.2,y0=0.85,y1=0.15,lwd=2,length=0.1)
    arrows(x0=1.6,y0=0.05,x1=0.4,lwd=2,length=0.1)
    arrows(x0=0.4,y0=0.1,x1=1.6,lwd=2,length=0.1)
    arrows(x0=1.8,y0=0.15,y1=0.85,lwd=2,length=0.1)
    arrows(x0=1.85,y0=0.85,y1=0.15,lwd=2,length=0.1)
    arrows(x0=1.6,y0=0.9,x1=0.4,lwd=2,length=0.1)
    arrows(x0=0.4,y0=0.95,x1=1.6,lwd=2,length=0.1)
    text(x=0.175,y=0.95,$(object.label[1]))
    text(x=1.825,y=0.95,$(object.label[2]))
    text(x=1.825,y=0.05,$(object.label[4]))
    text(x=0.175,y=0.05,$(object.label[3]))
    """
    R"""
    text(x=1,y=1,round($(object.rate[5]),signif),cex=0.8)
    """
    R"""
    text(x=1,y=0.85,round($(object.rate[6]),signif),cex=0.8)
    text(x=1.9,y=0.5,round($(object.rate[3]),signif),cex=0.8,srt=90)
    text(x=1.75,y=0.5,round($(object.rate[4]),signif),cex=0.8,srt=90)
    """
    R"""
    text(x=1,y=0,round($(object.rate[8]),signif),cex=0.8)
    text(x=1,y=0.15,round($(object.rate[7]),signif),cex=0.8)
    text(x=0.1,y=0.5,round($(object.rate[2]),signif),cex=0.8,srt=90)
    text(x=0.25,y=0.5,round($(object.rate[1]),signif),cex=0.8,srt=90)
    """
    return nothing
end
