import {Component, ElementRef, ViewChild} from '@angular/core';

@Component({
  selector: 'app-igv',
  templateUrl: './igv.component.html',
  styleUrls: ['./igv.component.scss']
})
export class IgvComponent {
  @ViewChild('igvDiv') igvDiv!: ElementRef;

  ngAfterViewInit(): void {
    setTimeout(() => {
      const options =
        {
          genome: "hg38",
          locus: "chr8:127,736,588-127,739,371",
          tracks: [
            {
              "name": "HG00103",
              "url": "https://s3.amazonaws.com/1000genomes/data/HG00103/alignment/HG00103.alt_bwamem_GRCh38DH.20150718.GBR.low_coverage.cram",
              "indexURL": "https://s3.amazonaws.com/1000genomes/data/HG00103/alignment/HG00103.alt_bwamem_GRCh38DH.20150718.GBR.low_coverage.cram.crai",
              "format": "cram"
            }
          ]
        };

      console.log(this.igvDiv);
      // @ts-ignore
      igv.createBrowser(this.igvDiv.nativeElement, options);
    });
  }
}
