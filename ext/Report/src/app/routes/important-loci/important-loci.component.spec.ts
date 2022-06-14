import { ComponentFixture, TestBed } from '@angular/core/testing';

import { ImportantLociComponent } from './important-loci.component';

describe('ImportantLociComponent', () => {
  let component: ImportantLociComponent;
  let fixture: ComponentFixture<ImportantLociComponent>;

  beforeEach(async () => {
    await TestBed.configureTestingModule({
      declarations: [ ImportantLociComponent ]
    })
    .compileComponents();
  });

  beforeEach(() => {
    fixture = TestBed.createComponent(ImportantLociComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
