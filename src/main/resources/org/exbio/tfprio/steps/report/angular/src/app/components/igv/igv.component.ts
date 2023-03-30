import {Component, ElementRef, ViewChild} from '@angular/core';
import {DataManagerService} from "../../services/data-manager.service";

@Component({
  selector: 'app-igv',
  templateUrl: './igv.component.html',
  styleUrls: ['./igv.component.scss']
})
export class IgvComponent {
  @ViewChild('igvDiv') igvDiv!: ElementRef;


  constructor(private dataService: DataManagerService) {
  }

  ngAfterViewInit(): void {
    setTimeout(() => {
      this.update();
    });
  }

  update(): void {
    this.igvDiv.nativeElement.innerHTML = '';

    const options =
      {
        genome: this.dataService.getConfigs()['InputConfigs']['genome'],
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

    // @ts-ignore
    igv.createBrowser(this.igvDiv.nativeElement, options);
  }
}
